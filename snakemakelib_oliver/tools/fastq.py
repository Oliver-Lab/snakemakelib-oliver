#!/usr/bin/env python
import gzip
import os
import re
import hashlib

import numpy as np

class FastqReader(object):
    """Returns a read-by-read fastQ parser analogous to file.readline() Based
    on: https://gist.github.com/xguse/1866279 """

    def __init__(self, filePath, headerSymbols=['@','+']):
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
            
    def __iter__(self):
        return self

    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            if isinstance(line, bytes):
                line = line.decode('utf-8')
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
        
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)

class FastqSummary(object):
    def __init__(self, fname):
        # Set up default counts
        self.fname = fname
        self.numReads = 0
        self._readLenArray = []
        self._header = []
        self._qualSet = set()
        self.notes = ''

        for read in FastqReader(fname):
            self.numReads += 1
            self.read = read
            self._readLenArray.append(len(read[1]))
            self._qualSet = self._qualSet.union(set(read[3]))
            self._header.append(read[0])

        self._parseHeader()
        self.readLen = self._checkReadLen()
        self.score = self._checkReadQual()
        self.md5 = self._md5sum()

    def _checkReadLen(self):
        rmin = np.min(self._readLenArray)
        rmax = np.max(self._readLenArray)
        if rmin == rmax:
            return rmin
        else:
            return '{min}:{max}'.format(rmin, rmax)

    def _checkReadQual(self):
        """
        https://en.wikipedia.org/wiki/FASTQ_format

          SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
          ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
          ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
          .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
          LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
          !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
          |                         |    |        |                              |                     |
         33                        59   64       73                            104                   126
          0........................26...31.......40                                
                                   -5....0........9.............................40 
                                         0........9.............................40 
                                            3.....9.............................40 
          0.2......................26...31........41                              

         S - Sanger        Phred+33,  raw reads typically (0, 40)
         X - Solexa        Solexa+64, raw reads typically (-5, 40)
         I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
         J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
             with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
             (Note: See discussion above).
         L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

        """

        values = np.array([ord(x) for x in self._qualSet])
        v = (np.min(values), np.max(values))

        # Remember python ranges don't include last number so you have to add 1

        # Sanger (0, 40) aka (33, 73)
        if v[0] in range(33, 73+1) and v[1] in range(33, 74+1):
            score = 'Phred+33 (Sanger)'

        # Illumina 1.8+ (0, 41) aka (33, 74)
        elif v[0] in range(0, 41+1) and v[1] in range(0, 41+1):
            score = 'Phred+33 (Illumina 1.8+)'
        
        # Illumina 1.5+ (3, 40) aka (67, 104)
        elif v[0] in range(67, 104+1) and v[1] in range(67, 104+1):
            score = 'Phred+64 (Illumina 1.5+)'

        # Illumina 1.3+ (0, 40) aka (64, 104)
        elif v[0] in range(64, 104+1) and v[1] in range(64, 104+1):
            score = 'Phred+64 (Illumina 1.3+)'
            
        # Solexa (-5, 40) aka (59, 104)
        elif v[0] in range(59, 104+1) and v[1] in range(59, 104+1):
            score = 'Phred+64 (Solexa)'

        else:
            self.note.append('ERROR: Quality scores appear to be outsitde of known ranges.')

        return score

    def _parseHeader(self):
        """
        https://en.wikipedia.org/wiki/FASTQ_format
        """ 
        self.instrument = ''
        self.runID = ''
        self.flowcellID = ''
        self.lane = ''
        self.pair = ''
        self.index = ''

        # grab the first read
        h = self._header[0]
        h = h.strip()

        # Pre CASAVA 1.8
        # @HWUSI-EAS100R:6:73:941:1973#0/1
        pattern = '(\@.*?):(\d+):\d+:\d+:\d+\#(.*)\/(\d)'
        if re.fullmatch(pattern, h):
            self.index = []
            for h in self._header:
                m = re.match(pattern, h)
                self.index.append(m.groups()[2])
            self.instrument = m.groups()[0]
            self.lane = m.groups()[1]
            self.pair = m.groups()[3]
            self.index = ':'.join([str(x) for x in set(self.index)])
            return

        # CASAVA 1.8
        # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
        pattern = '(\@.*?):(\d+):(.*?):(\d+):\d+:\d+:\d+ (\d):[Y,N]:\d+:(.*)'
        if re.fullmatch(pattern, h):
            self.index = []
            for h in self._header:
                m = re.match(pattern, h)
                self.index.append(m.groups()[5])
            self.instrument = m.groups()[0]
            self.runID = m.groups()[1]
            self.flowcellID = m.groups()[2]
            self.lane = m.groups()[3]
            self.pair = m.groups()[4]
            self.index = ':'.join([str(x) for x in set(self.index)])
            return

        # CASAVA 1.8 Truncated
        # @EAS139:136:FC706VJ:2:2104:15343:197393
        pattern = '(\@.*?):(\d+):(.*?):(\d+):\d+:\d+:\d+'
        if re.fullmatch(pattern, h):
            m = re.match(pattern, h)
            self.instrument = m.groups()[0]
            self.runID = m.groups()[1]
            self.flowcellID = m.groups()[2]
            self.lane = m.groups()[3]
            return

        # Probably SRA
        # @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
        #TODO: look at a bunch of SRAs and figure out pattern
        pattern = '(\@SRR.*?).*'
        if re.match(pattern, h):
            self.note.append("WARNINING: SRA header parsing is not implemented.")
            return

        self.note.append("ERROR: Header does not fit known patterns.")
        return

    def _md5sum(self):
        with open(self.fname, 'rb') as IN:
            md5 = hashlib.md5(IN.read()).hexdigest()
        return md5

    def getSummary(self):
        header = 'fileName,md5sum,numberReads,readLength,qualityScoreType,instrumentName,runID,flowcellID,laneNumber,matePair,index\n'
        res = [os.path.basename(self.fname), self.md5, self.numReads, self.readLen, 
               self.score, self.instrument, self.runID, self.flowcellID, self.lane, self.pair, self.index]
        return header + ','.join([str(x) for x in res])


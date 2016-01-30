# -*- snakemake -*-
from os.path import join
from snakemakelib.config import SNAKEMAKELIB_RULES_PATH
from snakemakelib_oliver.config import OLIVER_RULES_PATH, OLIVER_WORKFLOWS_PATH
from snakemake.utils import update_config

##############################
# Default configuration settings
##############################
qc_config = {
    'workflows.oliver.qc': {
        'aligner': 'bowtie2',
        'peakcallers': ['dfilter', 'macs2'],
        'trimadaptor': True,
        'bamfilter': True,
    },
    'settings': {
        'temp_rules': [],
    },
    'bio.ngs.settings': {
        'fastq_suffix': '.fastq.gz',
        'align_suffix': '.bam',
        'threads': 1,
    },
    'tmpdir': '/tmp'
}

update_config(qc_config, config)
config = qc_config

## General Include statements
include: join(SNAKEMAKELIB_RULES_PATH, 'settings.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'utils.rules')
include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs", "settings.rules")

## COUNT THE NUMBER OF READS
include: join(OLIVER_RULES_PATH, "count.rules")

config['bio.ngs.settings']['filter_suffix'] = ''
COUNT_TARGETS = generic_target_generator(tgt_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         target_suffix = ".fastq.count", 
                                         src_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         **config['bio.ngs.settings'])


## RUN FASTQC ON FASTQ
include: join(OLIVER_RULES_PATH, "qc", "fastqc.rules")

config['bio.ngs.settings']['filter_suffix'] = ''
FASTQC_TARGETS = generic_target_generator(tgt_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         target_suffix = "_fastqc_report.html", 
                                         src_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         **config['bio.ngs.settings'])


## RUN BASIC ALIGNMENT
aln_config = {'bio.ngs.align.bowtie2' : {
                 'cmd': 'module load bowtie/2-2.2.6; module load samtools/1.2; bowtie2',
                 'align': {
                     'threads' : config['bio.ngs.settings']['threads'],
                     'options' : '',
                 },
             },
           }

update_config(aln_config, config)
config = aln_config

include: join(SNAKEMAKELIB_RULES_PATH, "bio/ngs/align", "bowtie2.rules")
ruleorder: bowtie2_align_se > bowtie2_align_pe

config['bio.ngs.settings']['filter_suffix'] = ''
ALIGN_TARGETS = generic_target_generator(tgt_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         target_suffix = ".bam", 
                                         src_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         **config['bio.ngs.settings'])

## RUN FASTQC ON BAM
include: join(OLIVER_RULES_PATH, "qc", "fastqc.rules")

config['bio.ngs.settings']['filter_suffix'] = ''
BAM_FASTQC_TARGETS = generic_target_generator(tgt_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         target_suffix = "_bam_fastqc_report.html", 
                                         src_re = config['bio.ngs.settings']['sampleorg'].run_id_re, 
                                         **config['bio.ngs.settings'])


### RULES FOR RUNNING QC WORKFLOW
rule run_qc:
    input: COUNT_TARGETS + FASTQC_TARGETS + ALIGN_TARGETS + BAM_FASTQC_TARGETS

rule qc_targets:
    """List currently defined targets"""
    run:
        print("READ COUNT: ", COUNT_TARGETS)
        print("FASTQC FASTQ: ", FASTQC_TARGETS)
        print("ALIGNMENTS: ", ALIGN_TARGETS)
        print("FASTQC BAM: ", BAM_FASTQC_TARGETS)

# -*- snakemake -*-
import os
from glob import glob
from os.path import join, dirname, relpath, exists, basename
import shutil

import pandas as pd

from snakemake.utils import update_config, set_temporary_output, set_protected_output
from snakemake.workflow import workflow

from snakemake_rules import SNAKEMAKE_RULES_PATH

from snakemakelib.io import make_targets, IOTarget, IOAggregateTarget
from snakemakelib.sample.input import initialize_input, convert_samples_to_list
from snakemakelib.application import SampleApplication, PlatformUnitApplication

from snakemakelib_oliver.rules import OLIVER_RULES_PATH

# Import app functions. Make sure to change this for new workflows
from snakemakelib_oliver.workflows.dge.app import *

# Collect information about inputs; _samples stores a list of
# dictionaries, where each dictionary contains information on a sample
config["samples"] = convert_samples_to_list(config.get("samples", None))
config["_samples"] = initialize_input(src_re=config['sample.settings']['sample_organization'].run_id_re, 
                                      sampleinfo=config['sample.settings'].get('sampleinfo', None), 
                                      metadata=config['sample.settings'].get('metadata', None), 
                                      metadata_filter=config['sample.settings'].get('metadata_filter', None), 
                                      filter_suffix=config['sample.settings'].get("filter_suffix", ""), 
                                      sample_column_map=config['sample.settings'].get("sample_column_map", ""), 
                                      sample_filter=config.get("samples", None))


##############################
# Input Function
##############################
# Add custom input functions for handling data, typically needed to make sure
# merging is done correctly. (Required for samtools merge rule)
def _dge_samtools_merge_bam_input_fn(wildcards):
    tgt = IOTarget(config['sample.settings']['sample_organization'].run_id_re.file, suffix='.merge.bam')
    m = config['sample.settings']['sample_organization'].sample_re.match(wildcards.prefix)
    l = [x for x in _samples if x['SM'] == m.groupdict()['SM']]
    return make_targets(tgt_re=tgt, samples=l)

def _rg_fn(prefix):
    """ function to generate read group """
    m = run_id_re.match(prefix)
    if 'PU' in m.groupdict():
        return '--rg-id {SM}_{PU} --rg-sample {SM} --rg-platform-unit {PU}'.format(SM=m.groupdict()['SM'], PU=m.groupdict()['PU'])
    else:
        return '--rg-id {SM} --rg-sample {SM}'.format(SM=m.groupdict()['SM'])


##############################
# Default Configuration
##############################
# Set configuration for current workflow, setting can be overridden by the
# user.
dge_config = {
    'sample.settings': {
        'stranded': 'True',
        'ERCC': 'True',
    },
    'dge.workflow' : {
        'aligner' : 'tophat',
        'aggregate_output_dir': 'aggregated_results',
        'report_output_dir': 'report',
    },
    'settings' : {
        'temporary_rules' : [],
        'protected_rules' : [],
        'tmpdir': '/tmp',
    },
    'comp.settings': {
        'python2': {
            'activate_cmd': ''
        },
    },
    'bio.ngs.settings': {
        'fastq_suffix': '.fastq.gz',
        'align_suffix': '.bam',
        'read1_label': '',
        'read2_label': ''
    },
    'bio.ngs.tools.samtools' : {
        'rules': [],
        'merge' : {
            'inputfun' : _dge_samtools_merge_bam_input_fn,
        },
        'sort': {
            'options': '-T {prefix}.tmp'
        },
    },
    'bio.ngs.rnaseq.tuxedo' : {
        'rules' : [],
        'version2': True,
        'rg_fn': _rg_fn,
        'cufflinks': {
            'options': ''
         },
        'tophat': {
            'threads': 8,
            'options': '--library-type fr-firststrand -G {gtf}'.format(gtf=config['bio.ngs.settings']['annotation']['transcript_annot_gtf']),
         },
    },
    'bio.ngs.db.ucsc' : {
        'rules' : ['ucsc_wigToBigWig', 'ucsc_bedGraphToBigWig', 'ucsc_fetchChromSizes'],
    },
    'bio.ngs.rnaseq.htseq': {
        'rules': [],
    },
    'bio.ngs.rnaseq.deseq2': {
        'rules': ['deseq2_LET', 'deseq2', 'deseq2_sampleTable']
    },
}
update_config(dge_config, config)
config = dge_config

# On many HPC environments, the "module load" system is used to load various
# versions of applications. If you are on such a system this adds the module
# load prefix.
#TODO: Make this part of a separate config file, so the user has more control over which versions are used.
if config['settings']['module_load']:
    module_config = {
        'bio.ngs.tools.samtools': {
            'cmd': 'module load samtools/1.2; samtools',
        },
        'bio.ngs.rnaseq.htseq': {
            'cmd': 'module load htseq/0.6.1p1; htseq-count',
        },
        'bio.ngs.rnaseq.tuxedo': {
            'cufflinks': {
                'cmd': 'module load cufflinks; cufflinks',
            },
            'tophat': {
                'cmd': 'module load tophat/2.1.0; tophat',
            },
        },
    }
    update_config(module_config, config)
    config = module_config

##############################
# Include statements
##############################

# Set up conditional rule import
## Htseq
if config['sample.settings']['stranded'].lower() == 'forward':
    config['bio.ngs.rnaseq.htseq']['rules'] = ['htseq_count_forward', 'htseq_count_reverse', 'htseq_count_nonstranded_intergenic', 'htseq_count_copy_forward_to_genic']
    _run_stranded = True
    _run_nonstranded = False
elif (config['sample.settings']['stranded'].lower() == 'reverse') or (config['sample.settings']['stranded'] is True):
    config['bio.ngs.rnaseq.htseq']['rules'] = ['htseq_count_forward', 'htseq_count_reverse', 'htseq_count_nonstranded_intergenic', 'htseq_count_copy_reverse_to_genic']
    _run_stranded = True
    _run_nonstranded = False
elif (config['sample.settings']['stranded'].lower() == 'nonstranded') or (config['sample.settings']['stranded'] is False):
    config['bio.ngs.rnaseq.htseq']['rules'] = ['htseq_count_nonstranded', 'htseq_count_nonstranded_intergenic', 'htseq_count_copy_nonstranded_to_genic']
    _run_stranded = False
    _run_nonstranded = True
else:
    raise Exception("you must set sample.settings['stranded'] to one of {True, False, forward, reverse, nonstranded}. True assumes reverse stranded and False assumes nonstranded.")

## Aligner
ALIGNER = config['dge.workflow']['aligner'] 
if ALIGNER == 'tophat':
    config['bio.ngs.settings']['align_prefix'] = '.tophat2'
    config['bio.ngs.rnaseq.tuxedo']['rules'] = ['tuxedo_cufflinks_from_bam', 'tuxedo_tophat_se']
    config['bio.ngs.tools.samtools']['rules'] = ['samtools_sort', 'samtools_index']
elif ALIGNER == 'bowtie2':
    config['bio.ngs.settings']['align_prefix'] = '.bowtie2'
    config['bio.ngs.tools.samtools']['rules'] = ['samtools_sort', 'samtools_index', 'samtools_sam2bam']
else:
    raise Exception("you must set dge.workflow['aligner'] to one of {bowtie2, tophat}")

ALIGN_PREFIX = config['bio.ngs.settings']['align_prefix']
ALIGN_SUFFIX = config['bio.ngs.settings']['align_suffix'] 

# Import all of the required rules.
if workflow._workdir is None:
    raise Exception("no workdir set, set workdir before include statement!")

include: join(SNAKEMAKE_RULES_PATH, 'utils.rules')
include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/tools', 'samtools.rules')
include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/rnaseq', 'tuxedo.rules')

include: join(OLIVER_RULES_PATH, 'utils.rules')
include: join(OLIVER_RULES_PATH, 'align', 'bowtie2.rules')
include: join(OLIVER_RULES_PATH, 'tools', 'counting.rules')
include: join(OLIVER_RULES_PATH, 'rnaseq', 'htseq.rules')
include: join(OLIVER_RULES_PATH, 'rnaseq', 'deseq2.rules')

localrules: deseq2, agg_align, dge_report, link_tophat2, \
             htseq_count_copy_forward_to_genic, htseq_count_copy_reverse_to_genic, \
             htseq_count_copy_nonstranded_to_genic

##############################
# Workflow-specific rules
##############################    
rule link_tophat2:
    input: "{prefix}.tophat2/accepted_hits.bam"
    output: "{prefix}.tophat2.bam"
    params: ""
    run:
        odir = os.path.dirname(output[0])
        oname = os.path.basename(output[0])
        rp = './' + os.path.relpath(str(input), odir)
        shell("cd {odir}; ln -s {rp} {output}".format(odir=odir, rp=rp, output=oname))

rule unique_aln:
    input: "{prefix}.bam"
    output: "{prefix}.unique.bam"
    params: cmd = config['bio.ngs.tools.samtools']['cmd']
    shell: "{params.cmd} view -h {input} | egrep '^@|NH:i:1[[:space:]]' | samtools view - -b > {output}"
    
##############################
# Targets
##############################
# Set up the file name targets, should not needed changed.
_samples = config["_samples"]
if not config.get("samples", None)  is None:
    _samples = [s for s in _samples if s["SM"] in config["samples"]]

raw_run_re = config['sample.settings']['sample_organization'].raw_run_re
run_id_re = config['sample.settings']['sample_organization'].run_id_re
sample_re = config['sample.settings']['sample_organization'].sample_re

_rg_fn('output/DamIDseq_1/DamIDseq_1_run1')

##############################
# Applications
##############################
apps = {}
# Applications are objects store input, output, postprocessing information for
# different rules. You need an application for the main rules of a workflow.

## Run level
### Alignment
Align = PlatformUnitApplication(
    name=ALIGNER,
    iotargets={
        'bam': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ALIGN_SUFFIX), None),
        'log': (IOTarget(run_id_re.file, suffix='.align.log'), 
                IOAggregateTarget(os.path.join(config['dge.workflow']['aggregate_output_dir'], "{}_summary.csv".format(ALIGNER))))
        },
    units=_samples,
    run=True
)
if ALIGNER == 'bowtie2':
    Align.register_post_processing_hook('log')(dge_bowtie2_post_processing_hook)
apps['Align'] = Align

### Counting
Counts = PlatformUnitApplication(
    name='counting',
    iotargets={
        'fq_count': (IOTarget(run_id_re.file, suffix='.fastq.count'), None),
        'bam_count': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.sort.bam.count'), None),
        },
    units=_samples,
    run=False
)
apps['Counts'] = Counts

## Sample Level
### Genic Coverage Counts
HtseqStranded = SampleApplication(
    name="htseq",
    iotargets={
        'forward': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".unique" + ".forward.htseq.counts"),
                        IOAggregateTarget(os.path.join(config['dge.workflow']['aggregate_output_dir'], "Forward.htseq.counts"))),
        'reverse': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".unique" + ".reverse.htseq.counts"),
                        IOAggregateTarget(os.path.join(config['dge.workflow']['aggregate_output_dir'], "Reverse.htseq.counts"))),
    },
    units=_samples,
    run=_run_stranded
)
HtseqStranded.register_aggregate_post_processing_hook('forward')(dge_htseq_post_processing_hook)
HtseqStranded.register_aggregate_post_processing_hook('reverse')(dge_htseq_post_processing_hook)

HtseqNonstranded = SampleApplication(
    name="htseq",
    iotargets={
        'nonstranded': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".unique" + ".nonstranded.htseq.counts"),
                        IOAggregateTarget(os.path.join(config['dge.workflow']['aggregate_output_dir'], "Nonstranded.htseq.counts"))),
    },
    units=_samples,
    run=_run_nonstranded
)
HtseqNonstranded.register_aggregate_post_processing_hook('nonstranded')(dge_htseq_post_processing_hook)

if config['sample.settings']['stranded']:
    apps['Htseq'] = HtseqStranded
else:
    apps['Htseq'] = HtseqNonstranded

### Intergenic Coverage Counts
HtseqInter = SampleApplication(
    name="htseq_intergenic",
    iotargets={
        'intergenic': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".unique" + ".intergenic.nonstranded.htseq.counts"),
                        IOAggregateTarget(os.path.join(config['dge.workflow']['aggregate_output_dir'], "Intergenic.nonstranded.htseq.counts"))),
    },
    units=_samples,
    run=True
)
HtseqInter.register_aggregate_post_processing_hook('intergenic')(dge_htseq_post_processing_hook)
apps['HtseqInter'] = HtseqInter

## Run DESeq2
DESeq2_TARGETS = [os.path.join(config['dge.workflow']['aggregate_output_dir'], "deseq2_results.csv"), os.path.join(config['dge.workflow']['aggregate_output_dir'], "deseq2_meta.txt")]
apps['DESeq2'] = DESeq2_TARGETS

## Report targets
#REPORT_TARGETS = [join(config['atacseq.workflow']['report_output_dir'], "atacseq_all_rulegraph.png"), join(config['atacseq.workflow']['report_output_dir'], "atacseq_summary.html")]
#REPORT_TARGETS = [join(config['qc.workflow']['report_output_dir'], "index.html"), ]
#BIGWIG_TARGETS = [x.replace(".bed", ".bed.wig.bw") for x in Dfilter.targets['bed']] + [x.replace("_peaks.xls", "_treat_pileup.bdg.bw") for x in Macs2.targets['xls']] + [x.replace("_peaks.xls", "_control_lambda.bdg.bw") for x in Macs2.targets['xls']]

##############################
# Collection rules
##############################
# Rules to pull everthing together and run.

rule run_dge:
    input: Align.targets['bam'] + \
           HtseqStranded.targets['forward'] + HtseqStranded.targets['reverse'] + HtseqNonstranded.targets['nonstranded'] + HtseqInter.targets['intergenic'] +\
           DESeq2_TARGETS

rule dge_align:
    input: Align.targets['bam']

rule dge_count:
    input: HtseqStranded.targets['forward'] + HtseqStranded.targets['reverse'] + HtseqNonstranded.targets['nonstranded'] + HtseqInter.targets['intergenic']

rule dge_deseq2:
    input: DESeq2_TARGETS

rule debug:
    input: Align.targets['bam'] +\
           HtseqStranded.targets['forward'] + HtseqStranded.targets['reverse'] + HtseqInter.targets['intergenic'] 

rule agg_align:
    input: align = Align.targets['log']
    output: align = Align.aggregate_targets['log']
    run:
        aggregate_results(Align)

rule agg_htseq:
    input: htseq_forward = HtseqStranded.targets['forward'],
           htseq_reverse = HtseqStranded.targets['reverse'],
           htseq_nonstranded = HtseqNonstranded.targets['nonstranded'],
           htseq_inter = HtseqInter.targets['intergenic']
    output: htseq_forward = HtseqStranded.aggregate_targets['forward'],
            htseq_reverse = HtseqStranded.aggregate_targets['reverse'],
            htseq_nonstranded = HtseqNonstranded.aggregate_targets['nonstranded'],
            htseq_inter = HtseqInter.aggregate_targets['intergenic']
    run:
        aggregate_results(HtseqStranded)
        aggregate_results(HtseqNonstranded)
        aggregate_results(HtseqInter)

rule dge_report:
    input: align = Align.aggregate_targets['log'],
           htseq_forward = HtseqStranded.aggregate_targets['forward'],
           htseq_reverse = HtseqStranded.aggregate_targets['reverse'],
           htseq_nonstranded = HtseqNonstranded.aggregate_targets['nonstranded'],
           rulegraph = join(config['dge.workflow']['report_output_dir'], "static/images/dge_all_rulegraph.png")
    output: html = join(config['dge.workflow']['report_output_dir'], "index.html")
    run:
        qc_summary(config, input, output, **apps)
    
# Set temporary and protected outputs
set_temporary_output(*[workflow.get_rule(x) for x in config['settings']['temporary_rules']])
set_protected_output(*[workflow.get_rule(x) for x in config['settings']['protected_rules']])

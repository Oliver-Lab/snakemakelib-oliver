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
from snakemakelib_oliver.workflows.qc.app import *

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
def _qc_samtools_merge_bam_input_fn(wildcards):
    tgt = IOTarget(config['sample.settings']['sample_organization'].run_id_re.file, suffix='.merge.bam')
    m = config['sample.settings']['sample_organization'].sample_re.match(wildcards.prefix)
    l = [x for x in _samples if x['SM'] == m.groupdict()['SM']]
    return make_targets(tgt_re=tgt, samples=l)


##############################
# Default Configuration
##############################
qc_config = {
    'qc.workflow' : {
        'aligner' : 'bowtie2',
        'fastqc' : True,
        'trimadaptor': True,
        'fastqc_trimmed': True,
        'fastqc_bam': True,
        'fastqc_trimmed_bam': True,
        'align_untrimmed': True,
        'align_trimmed': True,
        'aggregate_output_dir': 'qc_aggregated_results',
        'report_output_dir': 'qc_report',
    },
    'settings' : {
        'temporary_rules' : [],
        'protected_rules' : [],
        'tmpdir': '/tmp',
    },
    'bio.ngs.settings': {
        'fastq_suffix': '.fastq.gz',
        'align_suffix': '.bam',
        'read1_label': '',
        'read2_label': ''
    },
    'bio.ngs.tools.samtools' : {
        'rules': ['samtools_sort', 'samtools_index'],
        'merge' : {
            'inputfun' : _qc_samtools_merge_bam_input_fn,
        },
    },
    'bio.ngs.qc.cutadapt' : {
        'rules' : ['cutadapt_cut_single_end'],
    },
    'bio.ngs.align.bowtie2' : {
        'rules' : ['bowtie2_align_se'],
        'align': {
            'threads': 4
        },
    },
    'bio.ngs.rnaseq.htseq': {
        'rules': ['htseq_count_forward', 'htseq_count_reverse', ],
    },
    'bio.ngs.db.ucsc' : {
        'rules' : ['ucsc_wigToBigWig', 'ucsc_bedGraphToBigWig', 'ucsc_fetchChromSizes'],
    },
}
update_config(qc_config, config)
config = qc_config

# Add module load if on a system with modules
if config['settings']['module_load']:
    module_config = {
        'bio.ngs.qc.fastqc': {
            'cmd': 'module load fastqc; fastqc ',
        },
        'bio.ngs.qc.cutadapt': {
            'cmd': 'module load cutadapt; cutadapt ',
        },
        'bio.ngs.align.bowtie2': {
            'cmd': 'module load bowtie/2-2.2.6; bowtie2',
        },
        'bio.ngs.tools.samtools': {
            'cmd': 'module load samtools/1.2; samtools',
        },
        'bio.ngs.rnaseq.htseq': {
            'cmd': 'module load htseq/0.6.1p1; htseq-count',
        },
    }
    update_config(module_config, config)
    config = module_config

##############################
# Include statements
##############################
# Import all of the required rules.
if workflow._workdir is None:
    raise Exception("no workdir set, set workdir before include statement!")

#include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/qc', 'fastqc.rules')
include: join(SNAKEMAKE_RULES_PATH, 'utils.rules')
include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/tools', 'samtools.rules')

include: join(OLIVER_RULES_PATH, 'qc', 'fastqc.rules')
include: join(OLIVER_RULES_PATH, 'align', 'bowtie2.rules')
include: join(OLIVER_RULES_PATH, 'rnaseq', 'htseq.rules')

if config['qc.workflow']['trimadaptor']:
    include: join(SNAKEMAKE_RULES_PATH, "bio/ngs/qc", "cutadapt.rules")

## Pick which alinger to use
ALIGNER = config['qc.workflow']['aligner'] 
if ALIGNER == 'tophat':
    config['bio.ngs.settings']['align_prefix'] = '.tophat2'
elif ALIGNER == 'bowtie2':
    config['bio.ngs.settings']['align_prefix'] = '.bowtie2'
    include: join(OLIVER_RULES_PATH, 'utils.rules')
    include: join(OLIVER_RULES_PATH, 'align', 'bowtie2.rules')
    include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/tools/_samtools', 'samtools_sam2bam.rule')
else:
    raise Exception("you must set qc.workflow['aligner'] to one of {bowtie2, tophat}")

ALIGN_PREFIX = config['bio.ngs.settings']['align_prefix']
ALIGN_SUFFIX = config['bio.ngs.settings']['align_suffix'] 

##############################
# Workflow-specific rules
##############################    


##############################
# Targets
##############################
_samples = config["_samples"]
if not config.get("samples", None)  is None:
    _samples = [s for s in _samples if s["SM"] in config["samples"]]

raw_run_re = config['sample.settings']['sample_organization'].raw_run_re
run_id_re = config['sample.settings']['sample_organization'].run_id_re
sample_re = config['sample.settings']['sample_organization'].sample_re

##############################
# Applications
##############################
apps = {}
# Run level targets
## FASTQC on FASTQ Files
Fastqc = PlatformUnitApplication(
    name="fastqc",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix="_fastqc/fastqc_report.html"), None),
        'Adapter Content': (IOTarget(run_id_re.file, suffix='_fastqc/Images/adapter_content.png'), None),
        'Duplication Levels': (IOTarget(run_id_re.file, suffix='_fastqc/Images/duplication_levels.png'), None),
        'Kmer Profiles': (IOTarget(run_id_re.file, suffix='_fastqc/Images/kmer_profiles.png'), None),
        'Per Base N Content': (IOTarget(run_id_re.file, suffix='_fastqc/Images/per_base_n_content.png'), None),
        'Per Base Quality': (IOTarget(run_id_re.file, suffix='_fastqc/Images/per_base_quality.png'), None),
        'Per Base Sequence Content': (IOTarget(run_id_re.file, suffix='_fastqc/Images/per_base_sequence_content.png'), None),
        'Per Sequence GC Content': (IOTarget(run_id_re.file, suffix='_fastqc/Images/per_sequence_gc_content.png'), None),
        'Per Sequence Quality': (IOTarget(run_id_re.file, suffix='_fastqc/Images/per_sequence_quality.png'), None),
        'Per Tile Quality': (IOTarget(run_id_re.file, suffix='_fastqc/Images/per_tile_quality.png'), None),
        'Sequence Length Distribution': (IOTarget(run_id_re.file, suffix='_fastqc/Images/sequence_length_distribution.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix='_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc']
)

Fastqc.register_post_processing_hook('summary')(qc_fastqc_summary_hook)
apps['Fastqc'] = Fastqc


## FASTQC on Trimmed FASTQ Files
Fastqc_trimmed = PlatformUnitApplication(
    name="fastqc_trimmed",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix=".trimmed_fastqc/fastqc_report.html"), None),
        'Adapter Content': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/adapter_content.png'), None),
        'Duplication Levels': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/duplication_levels.png'), None),
        'Kmer Profiles': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/kmer_profiles.png'), None),
        'Per Base N Content': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/per_base_n_content.png'), None),
        'Per Base Quality': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/per_base_quality.png'), None),
        'Per Base Sequence Content': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/per_base_sequence_content.png'), None),
        'Per Sequence GC Content': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/per_sequence_gc_content.png'), None),
        'Per Sequence Quality': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/per_sequence_quality.png'), None),
        'Per Tile Quality': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/per_tile_quality.png'), None),
        'Sequence Length Distribution': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/Images/sequence_length_distribution.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix='.trimmed_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_trimmed_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc_trimmed'] 
)

Fastqc_trimmed.register_post_processing_hook('summary')(qc_fastqc_summary_hook)
apps['Fastqc_trimmed'] = Fastqc_trimmed


## FASTQC on BAM Files
Fastqc_bam = PlatformUnitApplication(
    name="fastqc_bam",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".bam_fastqc/fastqc_report.html"), None),
        'Adapter Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/adapter_content.png'), None),
        'Duplication Levels': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/duplication_levels.png'), None),
        'Kmer Profiles': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/kmer_profiles.png'), None),
        'Per Base N Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/per_base_n_content.png'), None),
        'Per Base Quality': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/per_base_quality.png'), None),
        'Per Base Sequence Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/per_base_sequence_content.png'), None),
        'Per Sequence GC Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/per_sequence_gc_content.png'), None),
        'Per Sequence Quality': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/per_sequence_quality.png'), None),
        'Per Tile Quality': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/per_tile_quality.png'), None),
        'Sequence Length Distribution': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/Images/sequence_length_distribution.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.bam_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_bam_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc_bam'] 
)

Fastqc_bam.register_post_processing_hook('summary')(qc_fastqc_summary_hook)
apps['Fastqc_bam'] = Fastqc_bam


## FASTQC on Trimmed BAM Files
Fastqc_trimmed_bam = PlatformUnitApplication(
    name="fastqc_trimmed_bam",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".trimmed.bam_fastqc/fastqc_report.html"), None),
        'Adapter Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/adapter_content.png'), None),
        'Duplication Levels': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/duplication_levels.png'), None),
        'Kmer Profiles': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/kmer_profiles.png'), None),
        'Per Base N Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/per_base_n_content.png'), None),
        'Per Base Quality': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/per_base_quality.png'), None),
        'Per Base Sequence Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/per_base_sequence_content.png'), None),
        'Per Sequence GC Content': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/per_sequence_gc_content.png'), None),
        'Per Sequence Quality': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/per_sequence_quality.png'), None),
        'Per Tile Quality': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/per_tile_quality.png'), None),
        'Sequence Length Distribution': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/Images/sequence_length_distribution.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed.bam_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_trimmed_bam_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc_trimmed_bam'] 
)

Fastqc_trimmed_bam.register_post_processing_hook('summary')(qc_fastqc_summary_hook)
apps['Fastqc_trimmed_bam'] = Fastqc_trimmed_bam


## Cutadapt
Cutadapt = PlatformUnitApplication(
    name="cutadapt",
    iotargets={
        'cutadapt':(IOTarget(run_id_re.file, suffix=".cutadapt_metrics"),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "cutadapt.metrics")))},
    units=_samples,
    run=config['qc.workflow']['trimadaptor']
)

Cutadapt.register_post_processing_hook('cutadapt')(qc_cutadapt_post_processing_hook)
Cutadapt.register_plot('cutadapt')(qc_cutadapt_plot_metrics)
apps['Cutadapt'] = Cutadapt


# Alignments
#FIXME This is just a hack to make bowtie2 logs work.
logNames = {
    'bowtie2': '.bwt2',
    'tophat2': '.tophat2'
}
## Align raw reads
Align = PlatformUnitApplication(
    name=ALIGNER,
    iotargets={
        'bam': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ALIGN_SUFFIX), None),
        'log': (IOTarget(run_id_re.file, suffix="{}.log".format(logNames[ALIGNER])), IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "{}_summary.csv".format(ALIGNER))))
        },
    units=_samples,
    run=config['qc.workflow']['align_untrimmed']
)
if ALIGNER == 'bowtie2':
    Align.register_post_processing_hook('log')(qc_bowtie2_post_processing_hook)
apps['Align'] = Align


## Align trimmed reads
Align_trimmed = PlatformUnitApplication(
    name=ALIGNER + ' trimmed',
    iotargets={
        'bam': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + '.trimmed' + ALIGN_SUFFIX), None),
        'log': (IOTarget(run_id_re.file, suffix=".trimmed{}.log".format(logNames[ALIGNER])), IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "{}_trimmed_summary.csv".format(ALIGNER))))
        },
    units=_samples,
    run=config['qc.workflow']['align_trimmed']
)

Align_trimmed.register_post_processing_hook('log')(qc_bowtie2_post_processing_hook)
apps['Align_trimmed'] = Align_trimmed

# Coverage Counts
Htseq = SampleApplication(
    name="htseq",
    iotargets={
        'forward': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".forward.htseq.counts"),
                        IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "QC.Forward.htseq.counts"))),
        'reverse': (IOTarget(run_id_re.file, suffix=ALIGN_PREFIX + ".reverse.htseq.counts"),
                        IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "QC.Reverse.htseq.counts"))),
    },
    units=_samples,
)
Htseq.register_post_processing_hook('forward')(qc_htseq_post_processing_hook)
apps['Htseq'] = Htseq


# Report targets
REPORT_TARGETS = [join(config['qc.workflow']['report_output_dir'], "index.html"), ]
#BIGWIG_TARGETS = [x.replace(".bed", ".bed.wig.bw") for x in Dfilter.targets['bed']] + [x.replace("_peaks.xls", "_treat_pileup.bdg.bw") for x in Macs2.targets['xls']] + [x.replace("_peaks.xls", "_control_lambda.bdg.bw") for x in Macs2.targets['xls']]

##############################
# Collection rules
##############################

rule run_qc:
    input: Fastqc.targets['report'] + Fastqc_trimmed.targets['report'] + Fastqc_bam.targets['report'] + Fastqc_trimmed_bam.targets['report'] + \
           Cutadapt.targets['cutadapt'] + \
           Align.targets['bam'] + Align_trimmed.targets['bam'] + \
           Htseq.targets['forward'] + Htseq.targets['reverse'] 

rule run_align:
    input: Fastqc_bam.targets['report'] + Fastqc_trimmed_bam.targets['report']

rule qc_agg_fastqc:
    input: fastqc = Fastqc.targets['summary'],
           fastqc_trimmed = Fastqc_trimmed.targets['summary'],
           fastqc_bam = Fastqc_bam.targets['summary']
    output: fastqc = Fastqc.aggregate_targets['summary'],
            fastqc_trimmed = Fastqc_trimmed.aggregate_targets['summary'],
            fastqc_bam = Fastqc_bam.aggregate_targets['summary']
    run:
        aggregate_results(Fastqc)
        aggregate_results(Fastqc_trimmed)
        aggregate_results(Fastqc_bam)
        aggregate_results(Fastqc_trimmed_bam)

rule qc_agg_cutadapt:
    input: cutadapt = Cutadapt.targets['cutadapt']
    output: cutadapt = Cutadapt.aggregate_targets['cutadapt']
    run:
        aggregate_results(Cutadapt)

rule qc_agg_align:
    input: align = Align.targets['log'],
           align_trimmed = Align_trimmed.targets['log']
    output: align = Align.aggregate_targets['log'],
            align_trimmed = Align_trimmed.aggregate_targets['log']
    run:
        aggregate_results(Align)
        aggregate_results(Align_trimmed)

rule qc_agg_htseq:
    input: htseq_forward = Htseq.targets['forward'],
           htseq_reverse = Htseq.targets['reverse']
    output: htseq_forward = Htseq.aggregate_targets['forward'],
            htseq_reverse = Htseq.aggregate_targets['reverse']
    run:
        aggregate_results(Htseq)

rule qc_report:
    input: fastqc = Fastqc.aggregate_targets['summary'],
           fastqc_trimmed = Fastqc_trimmed.aggregate_targets['summary'],
           fastqc_bam = Fastqc_bam.aggregate_targets['summary'],
           cutadapt = Cutadapt.aggregate_targets['cutadapt'],
           align = Align.aggregate_targets['log'],
           align_trimmed = Align_trimmed.aggregate_targets['log'],
           htseq_forward = Htseq.aggregate_targets['forward'],
           htseq_reverse = Htseq.aggregate_targets['reverse'],
           rulegraph = join(config['qc.workflow']['report_output_dir'], "static/images/run_qc_rulegraph.png")
    output: html = join(config['qc.workflow']['report_output_dir'], "index.html")
    run:
        qc_summary(config, input, output, apps)
    
# Set temporary and protected outputs
set_temporary_output(*[workflow.get_rule(x) for x in config['settings']['temporary_rules']])
set_protected_output(*[workflow.get_rule(x) for x in config['settings']['protected_rules']])

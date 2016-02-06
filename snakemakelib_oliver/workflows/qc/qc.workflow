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
        'bamfilter' : True,
        'aggregate_output_dir': 'aggregated_results',
        'report_output_dir': 'report',
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
    'bio.ngs.tools.bamtools' : {
        'filter' : {
            'options' : {'mapQuality': ">=30",
                         'isProperPair': 'true'},
        },
    },
    'bio.ngs.qc.cutadapt' : {
        'rules' : ['cutadapt_cut_single_end'],
    },
    'bio.ngs.align.bowtie2' : {
        'rules' : ['bowtie2_align_se'],
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
            'cmd': 'module load fastqc' + config['bio.ngs.qc.fastqc']['cmd'],
        },
        'bio.ngs.qc.cutadapt': {
            'cmd': 'module load cutadapt' + config['bio.ngs.qc.cutadapt']['cmd'],
        },
        'bio.ngs.align.bowtie2': {
            'cmd': 'module load bowtie/2-2.2.6' + config['bio.ngs.align.bowtie2']['cmd'],
        },
        'bio.ngs.tools.samtools': {
            'cmd': 'module load samtools/1.2' + config['bio.ngs.tools.samtools']['cmd'],
        },
    }
    update_config(module_config, config)
    config = module_config


##############################
# Include statements
##############################
#include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/qc', 'fastqc.rules')
include: join(OLIVER_RULES_PATH, 'qc', 'fastqc.rules')
include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/align', 'bowtie2.rules')

if config['qc.workflow']['trimadaptor']:
    include: join(SNAKEMAKE_RULES_PATH, "bio/ngs/qc", "cutadapt.rules")

if workflow._workdir is None:
    raise Exception("no workdir set, or set after include of 'atacseq.workflow'; set workdir before include statement!")


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

## FASTQC on BAM Files
Fastqc_bam = PlatformUnitApplication(
    name="fastqc_bam",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix=".bam_fastqc/fastqc_report.html"), None),
        'Adapter Content': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/adapter_content.png'), None),
        'Duplication Levels': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/duplication_levels.png'), None),
        'Kmer Profiles': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/kmer_profiles.png'), None),
        'Per Base N Content': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/per_base_n_content.png'), None),
        'Per Base Quality': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/per_base_quality.png'), None),
        'Per Base Sequence Content': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/per_base_sequence_content.png'), None),
        'Per Sequence GC Content': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/per_sequence_gc_content.png'), None),
        'Per Sequence Quality': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/per_sequence_quality.png'), None),
        'Per Tile Quality': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/per_tile_quality.png'), None),
        'Sequence Length Distribution': (IOTarget(run_id_re.file, suffix='.bam_fastqc/Images/sequence_length_distribution.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix='.bam_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_bam_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc_bam'] 
)

Fastqc_bam.register_post_processing_hook('summary')(qc_fastqc_summary_hook)

## FASTQC on Trimmed BAM Files
Fastqc_trimmed_bam = PlatformUnitApplication(
    name="fastqc_trimmed_bam",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix=".trimmed.bam_fastqc/fastqc_report.html"), None),
        'Adapter Content': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/adapter_content.png'), None),
        'Duplication Levels': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/duplication_levels.png'), None),
        'Kmer Profiles': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/kmer_profiles.png'), None),
        'Per Base N Content': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/per_base_n_content.png'), None),
        'Per Base Quality': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/per_base_quality.png'), None),
        'Per Base Sequence Content': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/per_base_sequence_content.png'), None),
        'Per Sequence GC Content': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/per_sequence_gc_content.png'), None),
        'Per Sequence Quality': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/per_sequence_quality.png'), None),
        'Per Tile Quality': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/per_tile_quality.png'), None),
        'Sequence Length Distribution': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/Images/sequence_length_distribution.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix='.trimmed.bam_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_trimmed_bam_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc_trimmed_bam'] 
)

Fastqc_trimmed_bam.register_post_processing_hook('summary')(qc_fastqc_summary_hook)

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

# Alignments with bowtie2
Align = PlatformUnitApplication(
    name=config['qc.workflow']['aligner'],
    iotargets={
        'bam': (IOTarget(run_id_re.file, suffix='.bam'),
                None)},
    units=_samples,
    run=config['qc.workflow']['align_untrimmed']
)

Align_trimmed = PlatformUnitApplication(
    name=config['qc.workflow']['aligner'],
    iotargets={
        'bam': (IOTarget(run_id_re.file, suffix='.trimmed.bam'),
                None)},
    units=_samples,
    run=config['qc.workflow']['align_trimmed']
)

# Report targets
#REPORT_TARGETS = [join(config['atacseq.workflow']['report_output_dir'], "atacseq_all_rulegraph.png"), join(config['atacseq.workflow']['report_output_dir'], "atacseq_summary.html")]
REPORT_TARGETS = [join(config['qc.workflow']['report_output_dir'], "qc_summary.html"), ]

#BIGWIG_TARGETS = [x.replace(".bed", ".bed.wig.bw") for x in Dfilter.targets['bed']] + [x.replace("_peaks.xls", "_treat_pileup.bdg.bw") for x in Macs2.targets['xls']] + [x.replace("_peaks.xls", "_control_lambda.bdg.bw") for x in Macs2.targets['xls']]

##############################
# Collection rules
##############################

rule qc_all:
    input: Fastqc.targets['report'] + Fastqc_trimmed.targets['report'] + Fastqc_bam.targets['report'] + Fastqc_trimmed_bam.targets['report'] +  \
           Cutadapt.targets['cutadapt'] + \
           Align.targets['bam'] + Align_trimmed.targets['bam'] 

rule qc_aggregate_fastqc_results:
    input: fastqc = Fastqc.targets['summary'],
           fastqc_trimmed = Fastqc_trimmed.targets['summary'],
           fastqc_bam = Fastqc_bam.targets['summary'],
    output: fastqc = Fastqc.aggregate_targets['summary'],
            fastqc_trimmed = Fastqc_trimmed.aggregate_targets['summary'],
            fastqc_bam = Fastqc_bam.aggregate_targets['summary']
    run:
        aggregate_results(Fastqc)
        aggregate_results(Fastqc_trimmed)
        aggregate_results(Fastqc_bam)
        aggregate_results(Fastqc_trimmed_bam)

rule qc_aggregate_cutadapt_results:
    input: cutadapt = Cutadapt.targets['cutadapt']
    output: fastqc = Cutadapt.aggregate_targets['cutadapt']
    run:
        aggregate_results(Cutadapt)

rule qc_report:
    input: fastqc = Fastqc.aggregate_targets['summary']
    output: html = join(config['qc.workflow']['report_output_dir'], "qc_summary.html")
    run:
        qc_summary(config, input, output, Fastqc, Fastqc_trimmed, Fastqc_bam, Fastqc_trimmed_bam, Cutadapt, Align, Align_trimmed)
    
##############################
# Workflow-specific rules
##############################    


# Set temporary and protected outputs
#set_temporary_output(*[workflow.get_rule(x) for x in config['settings']['temporary_rules']])
#set_protected_output(*[workflow.get_rule(x) for x in config['settings']['protected_rules']])

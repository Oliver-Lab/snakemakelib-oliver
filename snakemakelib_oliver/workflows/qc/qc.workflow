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
        'trimadaptor' : True,
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
    },
    'bio.ngs.qc.fastqc': {
        'options': '-q --extract',
    },
    'bio.ngs.tools.bamtools' : {
        'filter' : {
            'options' : {'mapQuality': ">=30",
                         'isProperPair': 'true'},
        },
    },
    'bio.ngs.qc.cutadapt' : {
        'rules' : ['cutadapt_cut_paired_end'] # , 'cutadapt_qc_summary'],
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
            'cmd': 'module load fastqc; fastqc',
        },
    }
    update_config(module_config, config)
    config = module_config


##############################
# Include statements
##############################
include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/qc', 'fastqc.rules')

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
FASTQC = PlatformUnitApplication(
    name="fastqc",
    iotargets={
        'report': (IOTarget(run_id_re.file, suffix="_fastqc/fastqc_report.html"), None),
        'images': (IOTarget(run_id_re.file, suffix='_fastqc/Images/*.png'), None),
        'summary': (IOTarget(run_id_re.file, suffix='_fastqc/summary.txt'),
                    IOAggregateTarget(os.path.join(config['qc.workflow']['aggregate_output_dir'], "fastqc_summary.txt")))},
    units=_samples,
    run=config['qc.workflow']['fastqc']
)

FASTQC.register_post_processing_hook('summary')(qc_fastqc_summary_hook)
#FATQC.register_plot('fastqc')(atacseq_cutadapt_plot_metrics)


# Misc targets
#REPORT_TARGETS = [join(config['atacseq.workflow']['report_output_dir'], "atacseq_all_rulegraph.png"), join(config['atacseq.workflow']['report_output_dir'], "atacseq_summary.html")]
REPORT_TARGETS = [join(config['qc.workflow']['report_output_dir'], "qc_summary.html"), ]

#BIGWIG_TARGETS = [x.replace(".bed", ".bed.wig.bw") for x in Dfilter.targets['bed']] + [x.replace("_peaks.xls", "_treat_pileup.bdg.bw") for x in Macs2.targets['xls']] + [x.replace("_peaks.xls", "_control_lambda.bdg.bw") for x in Macs2.targets['xls']]


##############################
# Collection rules
##############################

rule qc_all:
    input: FASTQC.targets['report'] + REPORT_TARGETS


rule qc_fastqc:
    input: FASTQC.targets['report']


rule qc_aggregate_fastqc_results:
    input: fastqc = FASTQC.targets['summary']
    output: fastqc = FASTQC.aggregate_targets['summary']
    run:
        aggregate_results(FASTQC)

rule qc_report:
    input: fastqc = FASTQC.aggregate_targets['summary']
    output: html = join(config['qc.workflow']['report_output_dir'], "qc_summary.html")
    run:
        qc_summary(config, input, output, FASTQC)
    
##############################
# Workflow-specific rules
##############################    


# Set temporary and protected outputs
#set_temporary_output(*[workflow.get_rule(x) for x in config['settings']['temporary_rules']])
#set_protected_output(*[workflow.get_rule(x) for x in config['settings']['protected_rules']])

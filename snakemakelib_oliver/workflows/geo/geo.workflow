# -*- snakemake -*-
import os
from os.path import join
from glob import glob

from snakemake import shell
from snakemake.utils import update_config, set_temporary_output, set_protected_output

from snakemakelib.io import make_targets, IOTarget, IOAggregateTarget
from snakemakelib.application import PlatformUnitApplication
from snakemakelib.sample.input import initialize_input, convert_samples_to_list
from snakemakelib_oliver.rules import OLIVER_RULES_PATH

# Import app functions. Make sure to change this for new workflows
from snakemakelib_oliver.workflows.geo.app import *

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

def _geo_input(wildcards):
    tgt = IOTarget(config['sample.settings']['sample_organization'].run_id_re.file, suffix=config['bio.ngs.settings']['fastq_suffix'])
    prefix = config['geo.workflow']['geoDir'] +  wildcards.prefix
    m = config['geo.workflow']['geo_re'].match(prefix)
    return tgt.fmt.format(**m.groupdict())


def _geo_htseq_input(wildcards):
    #TODO: Need to fix the tophat2.unique part to work with the config system.
    tgt = IOTarget(config['sample.settings']['sample_organization'].run_id_re.file, suffix='.tophat2.unique.reverse.htseq.counts')
    prefix = config['geo.workflow']['geoDir'] +  wildcards.prefix
    m = config['geo.workflow']['geo_re'].match(prefix)
    return tgt.fmt.format(**m.groupdict())
    
##############################
# Default Configuration
##############################
# Set configuration for current workflow, setting can be overridden by the
# user.
geo_config = {
    'sample.settings': {
    },
    'geo.workflow' : {
        'geoDir': 'geo',
        'aggregate_output_dir': 'geo',
        'geo_re': IOTarget(join(config['geo.workflow']['geoDir'], "{SM}_{PU,[a-zA-Z0-9_]+}")),
    },
    'settings' : {
        'temporary_rules' : [],
        'protected_rules' : [],
        'tmpdir': '/tmp',
    }
}
update_config(geo_config, config)
config = geo_config


##############################
# Include statements
##############################
include: join(OLIVER_RULES_PATH, 'tools', 'fastqSummary.rules')

##############################
# Workflow-specific rules
##############################    
rule symlinkGeo:
    input: _geo_input
    output: config['geo.workflow']['geoDir'] + '{prefix}' + config['bio.ngs.settings']['fastq_suffix']
    run:
        from snakemakelib_oliver.utils import relative_symlink
        relative_symlink(str(input), str(output))

rule md5Counts:
    input: _geo_htseq_input
    output: config['geo.workflow']['geoDir'] + '{prefix}' + 'reverse.htseq.counts.md5sum'
    run:
        from snakemakelib_oliver.utils import file_md5sum
        md5 = file_md5sum(str(input))
        with open(output[0], 'w') as OUT:
            OUT.write('fileName,fileType,fileChecksum\n')
            OUT.write(','.join([str(input), 'text', str(md5)]))


localrules: symlinkGeo, md5Counts
##############################
# Targets
##############################
# Set up the file name targets, should not needed changed.
_samples = config["_samples"]
if not config.get("samples", None) is None:
    _samples = [s for s in _samples if s["SM"] in config["samples"]]

##############################
# Applications
##############################
apps = {}
# Applications are objects store input, output, postprocessing information for
# different rules. You need an application for the main rules of a workflow.

## Run level
### Alignment
Link = PlatformUnitApplication(
    name='Link',
    iotargets={'fastq': (IOTarget(config['geo.workflow']['geo_re'].file, suffix=config['bio.ngs.settings']['fastq_suffix']), None)},
    units=_samples,
    run=True
)
apps['Link'] = Link

RawFiles = PlatformUnitApplication(
    name='RawFiles',
    iotargets={'summary': (IOTarget(config['geo.workflow']['geo_re'].file, suffix='.fastq.summary'),
                           IOAggregateTarget(os.path.join(config['geo.workflow']['aggregate_output_dir'], "geo_summary.csv")))},
    units=_samples,
    run=True
)
RawFiles.register_post_processing_hook('summary')(geo_summary_table_post_processing_hook)
apps['RawFiles'] = RawFiles

CountFiles = PlatformUnitApplication(
    name='CountFiles',
    iotargets={'md5sum': (IOTarget(config['geo.workflow']['geo_re'].file, suffix='.reverse.htseq.counts.md5sum'),
                           IOAggregateTarget(os.path.join(config['geo.workflow']['aggregate_output_dir'], "geo_htseq_md5sum.csv")))},
    units=_samples,
    run=True
)
CountFiles.register_post_processing_hook('md5sum')(geo_summary_table_post_processing_hook)
apps['CountFiles'] = CountFiles

##############################
# Collection rules
##############################
# Rules to pull everthing together and run.
rule run_geo:
    input: Link.targets['fastq'] + RawFiles.targets['summary'] + CountFiles.targets['md5sum'] + \
           [RawFiles.aggregate_targets['summary']] + [CountFiles.aggregate_targets['md5sum']]

def aggregate_results(res):
    res.aggregate().save_aggregate_data()

rule agg_summary:
    input: summary = RawFiles.targets['summary']
    output: summary = RawFiles.aggregate_targets['summary']
    run:
        aggregate_results(RawFiles)

rule agg_md5sum:
    input: CountFiles.targets['md5sum']
    output: CountFiles.aggregate_targets['md5sum']
    run:
        aggregate_results(CountFiles)

# Set temporary and protected outputs
set_temporary_output(*[workflow.get_rule(x) for x in config['settings']['temporary_rules']])
set_protected_output(*[workflow.get_rule(x) for x in config['settings']['protected_rules']])

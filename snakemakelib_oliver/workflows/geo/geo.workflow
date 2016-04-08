# -*- snakemake -*-
import os
from glob import glob

from snakemake import shell

from snakemakelib.io import make_targets, IOTarget, IOAggregateTarget
from snakemakelib.sample.input import initialize_input, convert_samples_to_list
from snakemakelib_oliver.rules import OLIVER_RULES_PATH

# Collect information about inputs; _samples stores a list of
# dictionaries, where each dictionary contains information on a sample
config["samples"] = convert_samples_to_list(config.get("samples", None))
config["_samples"] = initialize_input(src_re = config['sample.settings']['sample_organization'].raw_run_re,
                                      sampleinfo = config['sample.settings'].get('sampleinfo', None),
                                      metadata = config['sample.settings'].get('metadata', None),
                                      metadata_filter = config['sample.settings'].get('metadata_filter', None),
                                      filter_suffix = config['sample.settings'].get("filter_suffix", ""),
                                      sample_column_map = config['sample.settings'].get("sample_column_map", ""),
                                      sample_filter = config.get("samples", None))

_samples = config["_samples"]
raw_run_re = config['sample.settings']['sample_organization'].raw_run_re
run_id_re = config['sample.settings']['sample_organization'].run_id_re
sample_re = config['sample.settings']['sample_organization'].sample_re

raw = IOTarget(raw_run_re.file, suffix=config['bio.ngs.settings']['fastq_suffix'])
run = IOTarget(run_id_re.file, suffix=config['bio.ngs.settings']['fastq_suffix'])

def _geo_input(x):
    

rule symlinkGeo:
    input: 
    output: 'geo/{prefix}' + config['bio.ngs.settings']['fastq_suffix']
    run:
        from snakemakelib_oliver.utils import relative_symlink
        relative_symlink(str(input), str(output))
        



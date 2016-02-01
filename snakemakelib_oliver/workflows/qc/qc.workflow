# -*- snakemake -*-
import os
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
config["_samples"] = initialize_input(src_re = config['settings']['sample_organization'].run_id_re,
                                      sampleinfo = config['settings'].get('sampleinfo', None),
                                      metadata = config['settings'].get('metadata', None),
                                      metadata_filter = config['settings'].get('metadata_filter', None),
                                      filter_suffix = config['settings'].get("filter_suffix", ""),
                                      sample_column_map = config['settings'].get("sample_column_map", ""),
                                      sample_filter = config.get("samples", None))

_samples = config["_samples"]


include: join(SNAKEMAKE_RULES_PATH, 'bio/ngs/qc', 'fastqc.rules')

fastqc_tgt = IOTarget(config['settings']['sample_organization'].run_id_re.file, suffix='_fastqc/fastqc_report.html')
rule run_qc:
    input: make_targets(fastqc_tgt, _samples)

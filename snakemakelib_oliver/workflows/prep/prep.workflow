# -*- snakemake -*-
import os
from glob import glob
from os.path import join, dirname, relpath, exists, basename
import shutil

from snakemakelib.io import make_targets, IOTarget, IOAggregateTarget
from snakemakelib.sample.input import initialize_input, convert_samples_to_list

# Collect information about inputs; _samples stores a list of
# dictionaries, where each dictionary contains information on a sample
config["samples"] = convert_samples_to_list(config.get("samples", None))
config["_samples"] = initialize_input(src_re = config['settings']['sample_organization'].raw_run_re,
                                      sampleinfo = config['settings'].get('sampleinfo', None),
                                      metadata = config['settings'].get('metadata', None),
                                      metadata_filter = config['settings'].get('metadata_filter', None),
                                      filter_suffix = config['settings'].get("filter_suffix", ""),
                                      sample_column_map = config['settings'].get("sample_column_map", ""),
                                      sample_filter = config.get("samples", None))

_samples = config["_samples"]

# Combine samples if needed
def combine_bins(raw, run):
    if not os.path.exists(run):
        os.makedirs(os.path.dirname(run), exist_ok=True)
        if len(raw) == 1:
            raw = raw[0]
            relpath = os.path.relpath(raw, os.path.dirname(run))
            odir = os.path.dirname(run)
            oname = os.path.basename(run)
            shell('cd {odir}; ln -fs {relpath} ./{output} && touch -h {output}'.format(odir=odir, relpath=relpath, output=oname))
        else:
            for i in raw:
                shell('cat {0} >> {1}'.format(i, run))

for s in _samples:
    raw = IOTarget(config['settings']['sample_organization'].raw_run_re.file, 
                   suffix=config['bio.ngs.settings']['fastq_suffix'])

    tgt_raw = make_targets(tgt_re=raw, samples=[s,])[0]
    raw_files = glob(tgt_raw)

    run = IOTarget(config['settings']['sample_organization'].run_id_re.file, 
                   suffix=config['bio.ngs.settings']['fastq_suffix'])
    run_files = make_targets(tgt_re=run, samples=[s,])[0]

    combine_bins(raw_files, run_files)

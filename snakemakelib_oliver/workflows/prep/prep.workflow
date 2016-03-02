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

def combine(input, output):
    # Only run if output does not exist
    if not os.path.exists(output):
        # Create output directory if needed
        os.makedirs(os.path.dirname(output), exist_ok=True)

        # Scan files and epxand any regular expressions
        globs = []
        for x in input:
            globs.extend(glob(x))

        # If only one file then symlink, otherwise concat.
        if len(globs) == 1:
            input = globs[0]
            relpath = os.path.relpath(input, os.path.dirname(output))
            odir = os.path.dirname(output)
            oname = os.path.basename(output)
            shell('cd {odir}; ln -fs {relpath} ./{output} && touch -h {output}'.format(odir=odir, relpath=relpath, output=oname))
        else:
            for i in globs:
                cmd = 'cat {0} >> {1}'.format(i, output)
                shell(cmd)
 

# Combine To Run Level
for r in make_targets(tgt_re=run, samples=_samples):
    m = run_id_re.match(r)
    l = [x for x in _samples if x['SM'] == m.groupdict()['SM']]
    combine(make_targets(raw, samples=l), r)

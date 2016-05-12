#!/bin/bash

# Intended to be called from the top-level dir of the repo, this installs the
# latest tagged git releases of snakemakelib-core and snakemakelib-rules, and
# installs snakemakelib-oliver.

set -e

SNAKEMAKELIBCORE=0.1-alpha.5
SNAKEMAKELIBRULES=0.1-beta.5

# install this package
python setup.py develop

pip install -e git+https://github.com/percyfal/snakemakelib-core.git@$SNAKEMAKELIBCORE#egg=snakemakelib-core
pip install -e git+https://github.com/percyfal/snakemake-rules.git@$SNAKEMAKELIBRULES#egg=snakemake-rules

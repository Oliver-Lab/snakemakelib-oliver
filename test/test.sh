#!/bin/bash
set -e
cd /opt/lcdb/examples
conda create -n snake --file requirements.txt --channel bioconda --channel r -y
source activate snake
snakemake

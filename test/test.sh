#!/bin/bash

cd /opt/lcdb/examples
conda env create --file snakemake_env.yaml --channel bioconda --channel r
source activate snake

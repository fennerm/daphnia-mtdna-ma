#!/usr/bin/env bash
snakemake --immediate-submit --debug --cluster 'snakefile_sbatch.py {dependencies}' -j 384 "$@"

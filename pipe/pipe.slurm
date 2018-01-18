#!/usr/bin/env bash
SCHEDULER="$HOME/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/sbatch_scheduler/sbatch_scheduler.py"
snakemake --cluster "$SCHEDULER" -j 499 --rerun-incomplete "$@"

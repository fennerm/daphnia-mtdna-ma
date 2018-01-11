#!/usr/bin/env bash
SCHEDULER="$HOME/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/slurm_scheduler.py"
snakemake --cluster "$SCHEDULER" --cluster-config "sbatch_config.yml" -j 499 \
    --rerun-incomplete "$@"

#!/usr/bin/env bash
SCHEDULER="$HOME/fmacrae/megadaph.private/megadaph/pipe/snakemake_slurm_utils/schedule_sbatch.py"
snakemake --keep-going --cluster "$SCHEDULER {dependencies}" -j 499 \
    --rerun-incomplete --notemp --immediate-submit --latency-wait 20 "$@"

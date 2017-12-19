"""Snakescript for locally realigning sequences around indels"""

from snakemake import shell

shell("gatk -T RealignerTargetCreator -R {snakemake.input.ref} -I "
      "{snakemake.input.bam} -o {snakemake.output.intervals}")

shell("gatk -T IndelRealigner -R {snakemake.input.ref} -I "
      "{snakemake.input.bam} -o {snakemake.output.bam} --targetIntervals "
      "{snakemake.output.intervals} --maxReadsForRealignment "
      "{snakemake.params.max_realignment_reads} --maxPositionalMoveAllowed "
      "{snakemake.params.max_move} --maxReadsForConsensuses "
      "{snakemake.params.max_cons_reads} --log_to_file {snakemake.log} "
      "--maxConsensuses {snakemake.params.max_cons})

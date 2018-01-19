"""Snakescript for extracting mitochondrial reads from the whole genome
alignments"""

import os
from tempfile import NamedTemporaryFile as tmp

from plumbum import local
from snakemake import shell

def get_nuclear_csome_names(bam):
    """Get the names of the mtdna contigs from the combined reference"""
    list_csomes = local["scripts/list_csomes"]
    names = list_csomes(bam)
    names = names.split("\n")
    names = names[2:]
    names = list(filter(None, names))
    names = ' '.join(names)
    return names

csomes = get_nuclear_csome_names(snakemake.input[0])

# Extract nuclear aligned reads
shell(' '.join(["samtools view -bh {snakemake.input}", csomes,
                "> {snakemake.output.bam}"]))
shell("samtools index {snakemake.output.bam}")

"""Snakescript for extracting mitochondrial reads from the whole genome
alignments"""

import os
from tempfile import NamedTemporaryFile as tmp

from plumbum.cmd import list_nonempty_csomes
from snakemake import shell


def get_mitochondrial_csome_names(bam):
    """Get the names of the mtdna contigs from the combined reference"""
    faids = list_nonempty_csomes(bam).split("\n")
    faids = ['"' + faid + '"' for faid in faids]
    mt_faids = ' '.join(faids[0:2])
    return mt_faids

csomes = get_mitochondrial_csome_names(snakemake.input[0])

# Names for .bam intermediary files
bam = tmp(suffix='.bam').name
bai = bam + '.bai'

# Extract mitochondrial aligned reads
shell("extract_csome {snakemake.input[0]} " + csomes +
      " | samtools sort - > " + bam)
shell("samtools index " + bam)

# Convert to .fastq
# SamToFastq throws an error if it finds unpaired mates. We can just ignore it
# since we're not interested in them.
shell("picard SamToFastq I=" + bam + " F={snakemake.output.fwd_reads} "
      "F2={snakemake.output.rev_reads} "
      "UNPAIRED_FASTQ={snakemake.output.unpaired_reads} "
      "VALIDATION_STRINGENCY=SILENT")

# Clean up
os.remove(bai)

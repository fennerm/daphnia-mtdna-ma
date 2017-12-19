"""Snakescript for extracting mitochondrial reads from the whole genome
alignments"""

import os
from tempfile import NamedTemporaryFile as tmp

from Bio.SeqIO import parse
from snakemake import shell


def next_seq_id(fasta_iterator):
    """Get the next sequence ID from a seqIO fasta iterator object"""
    return str(next(fasta_iterator).id)


def get_mitochondrial_csome_names(combined_reference):
    """Get the names of the mtdna contigs from the combined reference"""
    seq_iter = parse(combined_reference, "fasta")
    faids = [next_seq_id(seq_iter), next_seq_id(seq_iter)]
    # The names need to be quoted to avoid bash expansion
    faids = ["'" + faid + "'" for faid in faids]
    faid_str = ' '.join(faids)
    return faid_str

csomes = get_mitochondrial_csome_names(snakemake.input.ref)

# Names for .bam intermediary files
bam = tmp(suffix='.bam').name
bai = bam + '.bai'

# Extract mitochondrial aligned reads
shell("samtools view -bh {snakemake.input.bam} " + csomes + " > " + bam)
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

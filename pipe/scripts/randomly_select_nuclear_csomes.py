#!/usr/bin/env python
import random

from plumbum import local
from plumbum.cmd import list_nonempty_csomes


def partition_by_species(samples):
    """Split the list of input .bam files by their species"""
    pulex = []
    magna = []
    for sample in samples:
        if sample.name.startswith(("F", "G", "I")):
            magna.append(sample)
        elif sample.name.startswith(("L", "TCO")):
            pulex.append(sample)
    return pulex, magna


samples = [local.path(sample) for sample in snakemake.input]
species = ("pulex", "magna")
species_split = partition_by_species(samples)
for spp, subsplit in zip(species, species_split):
    # Get all nuclear csomes with at least one alignment
    # The first two chromosomes of the .bam files are mitochondrial
    csomes = [
        list_nonempty_csomes(sample).split("\n")[2:] for sample in subsplit
    ]

    # Get csomes which have at least one alignment in all samples
    shared_csomes = set(csomes[0]).intersection(*csomes)

    # Randomly select some chromosomes
    if int(snakemake.params.n) < len(shared_csomes):
        rand_csomes = random.sample(shared_csomes, snakemake.params.n)
    else:
        # If number of shared csomes is less than n, just return them all
        # (This should only happen in testing)
        rand_csomes = shared_csomes

    # Convert to string and write to output file
    rand_csomes = '\n'.join(rand_csomes)
    outfile = snakemake.output[spp]
    with open(outfile, "w") as out:
        out.write(rand_csomes)

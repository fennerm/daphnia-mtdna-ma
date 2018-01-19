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
    # Get all csomes with at least one alignment
    csomes = [list_nonempty_csomes(sample).split("\n") for sample in subsplit]

    # Get csomes which have at least one alignment in all samples
    shared_csomes = set(csomes[0]).intersection(*csomes)

    # Randomly select some chromosomes
    rand_csomes = random.sample(shared_csomes, snakemake.params.n)

    # Convert to string and write to output file
    shared_csomes = '\n'.join(shared_csomes)
    outfile = snakemake.output[spp]
    with open(outfile) as out:
        out.write(shared_csomes)

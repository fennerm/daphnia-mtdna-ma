#!/usr/bin/env python

"""NOVOPlasty assembly.

Usage:
	novoplasty_assemble.py -s <reference_fasta> -l <read_length> <isolates>...
	novoplasty_assemble.py (-h | --help)

Options:
	-h --help			Show this screen.
	-s --seed			Reference sequence 
	-l --read_length	Integer; length of untrimmed reads
"""

## Produce NOVOPlasty mitochondrial assemblies from Megadaph data. 
## One assembly is produced per isolate
from __future__ import print_function                                           
from pathlib import Path
import pdb
import os, sys, contextlib
from docopt import docopt
from fen_util import (gzip, gunzip, mkdir, mkdirs, paths_to_string, 
	working_directory, concat)

# NOVOPLASTY Parameters
INSERT_SIZE = 320
INSERT_SIZE_AUT = 'yes'
TYPE = 'mito'
GENOME_RANGE = '15000-20000'
K_MER = 80
INSERT_RANGE = 1.6
INSERT_RANGE_STRICT = 1.1
SINGLE_PAIRED = 'PE'
EXTENDED_LOG = 0
SAVE_ASSEMBLED_READS = 'no'

# Bowtie2 parameters
MAX_INSERT_SIZE = 10000
THREADS = 40
PRESET = "--very-sensitive-local"
FLAGS = ""

# Write a NOVOPlasty config file for the given parameters
def create_config_file(genotype, seed, read_length, reads, out_dir):
	config_path = Path(out_dir, 'config.txt')

	with open(str(config_path), 'w') as f:
		print("Project name         = " + genotype, file = f)
		print("Insert size          = " + str(INSERT_SIZE), file=f)
		print("Insert size aut      = " + INSERT_SIZE_AUT, file=f)
		print("Read Length          = " + read_length, file=f)
		print("Type                 = " + TYPE, file=f)
		print("Genome Range         = " + GENOME_RANGE, file=f)
		print("K-mer                = " + str(K_MER), file=f)
		print("Insert Range         = " + str(INSERT_RANGE), file=f)
		print("Insert Range strict  = " + str(INSERT_RANGE_STRICT), file=f)
		print("Single/Paired        = " + SINGLE_PAIRED, file=f)
		print("Max memory           =", file=f)
		print("Extended log         = " + str(EXTENDED_LOG), file=f)
		print("Save assembled reads = " + SAVE_ASSEMBLED_READS, file=f)
		print("Combined reads       =", file=f)
		print("Forward reads        = " + str(reads[0]), file=f)
		print("Reverse reads        = " + str(reads[1]), file=f)
		print("Seed Input           = " + seed, file=f)
	
	return config_path

def novoplasty(genotype, genotype_dir, reads, read_length, seed):
	config_path = create_config_file(genotype, seed, read_length, reads, \
		genotype_dir)

	with working_directory(str(genotype_dir)):
		os.system("NOVOPlasty2.6.2.pl -c " + str(config_path))

def bowtie2(

if __name__ == "__main__":

	# Parse Args
	arg = docopt(__doc__)
	print(arg)
	isolates = arg['<isolates>']
	read_length = arg['<read_length>']
	seed = arg['<reference_fasta>']
	
	# Initialize directories
	root_dir = Path.cwd()
	cons_dir = Path(root_dir, 'create_consensus')
	mkdir(cons_dir)
	isolate_dirs = mkdirs(isolates, cons_dir)

	for i in range(len(isolates)):
		
		controls = sorted(Path(root_dir).glob(isolates[i] + '*SC_*.gz'))

		novoplasty(isolates[i], isolate_dirs[i], zipped, read_length, seed)
		#samples = samples + sorted(Path(root_dir).glob(isolates[i] + 
		#	'*EC*.gz'))
		#print(samples)

		## Unzip
		#for s in samples:
		#	gunzip(s, isolate_dirs[i])

		#zipped = []

		## Concatenate
		#for n in [1, 2]:
		#	reads = sorted(isolate_dirs[i].glob('*_'+ str(n) + '*'))
		#	concat_name = isolates[i] + '_' + str(n) + '.fastq'
		#	concat_out = Path(isolate_dirs[i], concat_name)
		#	concat(reads, concat_out)

		#	# Gzip the concatenated file
		#	gzip_out = gzip(concat_out)
		#	relative = gzip_out.relative_to(isolate_dirs[i])
		#	zipped.append(relative)

		#	# Remove unzipped
		#	for r in reads:
		#		Path.unlink(r)


#!/usr/bin/env python
""" Pipeline for the mtDNA analysis portion of Megadaph

Directory locations are hardlinked based on the species and isolate ID so can
be run from any directory.

Not particularly proud of the implmentation of this one - was written during
my undergraduate thesis and not much thought was put into extensibility or
maintenance. Would not recommend attempting to use any of this code for other
project but it did the job here.

Usage:
    daphnia_mt_ma_pipeline.py -s <species> -i <isolate>
    dapnia_mt_ma_pipeline.py (-h | --help)

Options:
    -h --help  Show this screen.
    -s --species  'pulex' or 'magna'
    -i --isolate  One of ['LIN', 'TCO', 'I', 'G', 'F']

"""
import os, shutil, sys, subprocess, multiprocessing, math
from docopt import docopt
from joblib import Parallel, delayed
from glob import glob
from psutil import virtual_memory


def init_dirs(directories):
    """Given a dictionary of directories.
    Create any that do not already exist."""
	for d, path in directories.items():
		if not os.path.isdir(path):
			os.makedirs(path)

def pair_files(directory):
    """Given a directory with paired fastq or fastq.gz files. Return a list of
    tuples of paired files."""
	p1 = (glob(os.path.join(directory, '*_1.*fq')) +
			glob(os.path.join(directory, '*_1.*fastq')))

	p2 = (glob(os.path.join(directory, '*_2.*fq')) +
			glob(os.path.join(directory, '*_2.*fastq')))

	p1_base = []
	for p in p1:
		p1_base.append(os.path.basename(p))
	p1_base = sorted(p1_base)

	p2_base = []
	for p in p2:
		p2_base.append(os.path.basename(p))
	p2_base = sorted(p2_base)

	return zip(p1_base, p2_base)

def sample_name_with_orient(f):
    """Return the sample name of a file, including the paired end orientation.
    """
	return f.split('.')[0]

def sample_name_root(f):
    """Return the sample name of a file, excluding the paired end
    orientation."""
	return sample_name_with_orient(f)[:-2]

def clear_dirs(dir_list):
    """ Given a directory. Delete all regular files in directory."""
	# Check that directory is within the pipeline tree, to avoid accidentally
    # deleting important irreplaceable files.
	for d in dir_list:
		if d not in dirs.values():
			print('Exiting. Pipeline attempted to delete contents of '+d)
			sys.exit()

		# Check that files in root directory are not being targeted for
        # deletion.
		elif d == dirs['root']:
			print('Exiting. Pipeline attempted to delete raw sequence files')
			sys.exit()

		else:
			for f in os.listdir(d):
				if not os.path.isdir(os.path.join(d, f)):
					if os.path.islink(os.path.join(d, f)):
						os.unlink(os.path.join(d, f))
					else:
						os.remove(os.path.join(d, f))

def run_command(command):
    """Run a bash command in python, and return a tuple of (STDOUT, STDERR).
    Command input is a list.
    Haven't really figured out how to get piping to work properly, so if piping
    is necessary, use os.system() instead."""
	ps = subprocess.Popen(command, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
	out, err = ps.communicate()
	print(out)
	print(err)

	return out, err

def determine_threads(files):
    """Given a list of files, return the number of the threads
    Not guaranteed to prevent memory errors since it assumes that tools
    required RAM approximately equal to input file size. If a tool tries to
    allocate a multiple of input file size, system might crash."""
	files_memory = 0

	# Loop through files and determine their cumulative file size in bytes.
	for f in files:
		files_memory = files_memory + os.stat(f).st_size

	# Average file size
	avg_size = files_memory/(len(files))

	# Determine free memory
	mem_avail = virtual_memory().available

	# Number of threads which will use up large amount of available memory (but not all).
	p = math.ceil(((mem_avail) / (10 * avg_size)))

	# If all files could easily be held in RAM, use a thread for each file.
	if p > len(files):
		p = len(files)

	return p


def align(wd, odir, ix, logdir):
    """Align reads in working directory (wd) to reference index (ix) and write
    logs to logdir."""

	pairs = pair_files(wd)
	print(pairs)
	max_ins = 100000

	for p in pairs:

		out_sam = os.path.join(odir, sample_name_root(p[0])+'.sam')

		p0 = os.path.join(wd, p[0])
		p1 = os.path.join(wd, p[1])
		# Run Bowtie2
		bowtie_out, bowtie_err = run_command(['bowtie2', bowtie2_preset, '-x',
			ix, '-X', max_ins, '-p', threads, '-1', p0, '-2', p1,
			'-S', out_sam])

		## Write logfiles

		bowtie_log = os.path.join(logdir, p[0].replace(".fq", ".log"))

		with open(bowtie_log, 'a') as out_file:
			out_file.write('BOWTIE STANDARD OUT \n')
			if bowtie_out is not None:
				out_file.write(bowtie_out)

			out_file.write('\nBOWTIE STANDARD ERROR\n')
			if bowtie_err is not None:
				out_file.write(bowtie_err)

def add_read_groups(f, odir):
    """Add read group information to a bam file."""
	f_name = os.path.basename(f)
	_id = sample_name_with_orient(f)
	id_root = sample_name_root(f)
	out_bam = os.path.join(odir, f_name.replace(".bam",".grouped.bam"))
	run_command(['picard', 'AddOrReplaceReadGroups', 'I='+f, 'O='+out_bam,
		'RGLB='+id_root, 'RGPL=illumina', 'RGPU=unit1', 'RGSM='+_id])
	run_command(["samtools", "index", out_bam])
	os.remove(f)
	os.remove(f+'.bai')

def local_realign(f, odir, ref):
    """Locally realign a bam file with GATK"""
	f_name = os.path.basename(f)
	intervals = os.path.join(odir, f_name.replace('.bam', '.intervals'))
	out_bam = os.path.join(odir, f_name.replace('.bam', '.realign.bam'))
	run_command(['gatk', '-T', 'RealignerTargetCreator', '-R', ref, '-I', f,
        '-o', intervals])
	run_command(['gatk', '-T', 'IndelRealigner', '-R', ref, '-I', f,
		'--targetIntervals', intervals, '-o', out_bam])
	os.remove(f)
	os.remove(f+'.bai')

def mark_dups(f, odir, logdir):
    """Mark duplicates in a bam file with picard."""
	f_name = os.path.basename(f)
	out_bam = os.path.join(odir, f_name.replace('.bam', '.mrkdup.bam'))
	logfile = os.path.join(logdir, f_name.replace('.bam', '.txt'))
	run_command(['picard', 'MarkDuplicates', 'I='+f, 'O='+out_bam,
		'REMOVE_DUPLICATES='+removeDups, 'METRICS_FILE='+logfile])
	run_command(['samtools', 'index', out_bam])
	os.remove(f)
	os.remove(f.replace('.bam', '.bai'))

def to_bam_and_sort(wd):
	"""Convert sam to bam"""
	print("Converting SAM to BAM")
	for f in os.listdir(wd):
		if f.endswith('.sam'):
			in_sam = os.path.join(wd, f)
			out_bam = os.path.join(wd, f.replace('.sam', '.bam'))
			run_command(['samtools', 'view', '-bhT', mtdna_app_ref, '-@',
                threads, '-o', out_bam, in_sam])
			os.remove(in_sam)

	## Sort bam files
	print("Sorting BAM files")
	for f in os.listdir(wd):
		if f.endswith('.bam'):
			in_bam = os.path.join(wd, f)
			out_bam = os.path.join(wd, f.replace('.bam', '.sorted.bam'))
			run_command(['samtools', 'sort', '-o', out_bam, '-@', threads,
                '-m', '4G', in_bam])
			run_command(['samtools', 'index', out_bam])
			os.remove(in_bam)

def viterbi_realignment(f_path, wd, ref):
	out_bam = os.path.join(wd, f_path.replace(".bam", ".viterbi.bam"))
	run_command(["lofreq", "viterbi", "-f", ref, "-o", out_bam, f_path])
	os.remove(f_path)
	os.remove(f_path+".bai")

def post_map(wd, mrkduplogs, ref):
    """Run set of processes which follow alignment. Includes local realignment,
    duplicate removal, sorting etc.
    Args:
        cwd = directory with files to be processed.
        mrkduplogs = directory for markduplicates logfiles
        ref = reference sequence corresponding to bam files"""

	## Add read group information (required by GATK)
	# Not using the group fields correctly, but they're not actually used
    # anyway.

	print("Adding read group information")
	for f in os.listdir(wd):
		if f.endswith('.bam'):
			f_path = os.path.join(wd, f)
			add_read_groups(f_path, wd)

	print("Viterbi realignment")
	for f in os.listdir(wd):
		if f.endswith('.bam'):
			f_path = os.path.join(wd, f)
			viterbi_realignment(f_path, wd, ref)

	print("Sorting BAM files")
	for f in os.listdir(wd):
		if f.endswith('.bam'):
			in_bam = os.path.join(wd, f)
			out_bam = os.path.join(wd, f.replace('.bam', '.sorted.bam'))
			run_command(['samtools', 'sort', '-o', out_bam, '-@', threads,
                '-m', '4G', in_bam])
			run_command(['samtools', 'index', out_bam])
			os.remove(in_bam)

	## Local realign
	print("Local Realignment")
	for f in os.listdir(wd):
		if f.endswith('.bam'):
			f_path = os.path.join(wd, f)
			local_realign(f_path, wd, ref)

    ## Remove duplicates
	print("Removing duplicates")
	for f in os.listdir(wd):
		if f.endswith('.bam'):
			f_path = os.path.join(wd, f)
			mark_dups(f_path, wd, mrkduplogs)


# Produce readcounts showing allele frequencies at every locus in ref
def counts(f, odir, ref):
	f_name = os.path.basename(f)
	out_counts = os.path.join(odir, f_name.replace(".bam", ".counts"))

	# Write counts
	os.system("bam-readcount -w 1 -f "+ref+" "+f+" > "+out_counts)

# Convert all bam files in wd to paired fastq files and place in odir.
def bam_to_fastq(f, odir):
	f_name = os.path.basename(f)
	reads_basename = sample_name_with_orient(f_name)
	p1=os.path.join(odir, reads_basename+'_1.fq')
	p2=os.path.join(odir, reads_basename+'_2.fq')
	run_command(['picard', 'SamToFastq', 'I='+f, 'F='+p1, 'F2='+p2])

# Call variants for all bam files in wd using lofreq.
def call_variants(wd, ref, odir):
	## Assign Indel Qualities
	print("Calling variants")
	for f in os.listdir(wd):
		if f.endswith(".bam"):
			in_bam = os.path.join(wd, f)
			out_bam = os.path.join(odir, f.replace(".bam", ".indelqual.bam"))
			run_command(['lofreq', 'indelqual', '--dindel', '-f', ref, '-o', out_bam, in_bam])
			run_command(["samtools", "index", out_bam])

	## Call Variants
	for f in os.listdir(odir):
		if f.endswith(".bam"):
			in_bam = os.path.join(odir, f)
			out_vcf = os.path.join(odir, f.replace(".bam", ".vcf"))
			run_command(["lofreq", "call", "-l", region, "-f", ref, "-o",
				out_vcf, "--call-indels", in_bam])

	cwd = os.getcwd()
	os.chdir(odir)

	vcfs = list()
	for f in os.listdir(odir):
		if f.endswith('.vcf'):
			vcf_path = os.path.join(odir, f)
			vcfs.append(vcf_path)

	## Find shared variants

	os.system("multi_vcf_subtract.py "+' '.join(vcfs))

	os.chdir(cwd)


def run_pipeline(start):
	if start == dirs['root']:
		run_qc()
	elif start == dirs['trim']:
		run_comp_map()
	elif start == dirs['comp_map']:
		run_map_to_og_and_rot(True)
	elif start == dirs['mtdna_comp_map']:
		run_map_to_og_and_rot(False)
	elif start == dirs['map_rot_og']:
		run_call_variants()

def run_qc():

	clear_dirs([dirs['trim'], dirs['skewerlog'], dirs['trim_fastqc']])

	# Unzip .gz files
	for f in os.listdir(dirs['root']):
		if f.endswith(".gz"):
			gz_in = os.path.join(dirs['root'], f)
			fq_out = os.path.join(dirs['root'], f.replace(".gz", ""))
			os.system('gunzip -c '+gz_in+' > '+fq_out)

	pairs = pair_files(dirs['root'])
	for p in pairs:
		p0 = os.path.join(dirs['root'], p[0])
		p1 = os.path.join(dirs['root'], p[1])
		run_command(["fastq-sort.pl", p0, p1])
		p0_out = os.path.join(dirs['trim'], p[0].replace(".fq", ".sort.fq"))
		p1_out = os.path.join(dirs['trim'], p[1].replace(".fq", ".sort.fq"))
		shutil.move(p0.replace(".fq", ".sort.fq"), p0_out)
		shutil.move(p1.replace(".fq", ".sort.fq"), p1_out)
		os.remove(os.path.join(dirs['root'], p[0].replace(".fq", ".singletons.fq")))
		os.remove(os.path.join(dirs['root'], p[1].replace(".fq", ".singletons.fq")))

	## Reencode quality scores in pulex dataset
	if spp == "pulex":
		for f in os.listdir(dirs['trim']):
			if f.endswith(".fq"):
				in_fq = os.path.join(dirs['trim'], f)
				out_fq = os.path.join(dirs['trim'], f.replace(".fq", ".rePhred.fq"))
				os.system('seqtk1.2 seq -Q 64 -V ' +  in_fq + '> ' + out_fq)
				os.remove(in_fq)

	# We rename the magna files with "rePhred" for convenience
	elif spp == "magna":
		for f in os.listdir(dirs['trim']):
			if f.endswith(".fq"):
				in_fq = os.path.join(dirs['trim'], f)
				out_fq = os.path.join(dirs['trim'], f.replace(".fq", ".rePhred.fq"))
				shutil.move(in_fq, out_fq)

	pairs = pair_files(dirs['trim'])

	## Skewer
	for p in pairs:
		p0_in = os.path.join(dirs['trim'], p[0])
		p1_in = os.path.join(dirs['trim'], p[1])
		if spp == "pulex":
			run_command(['skewer', '-t', threads, '-q', qual_threshold, '-u', p0_in, p1_in])
		elif spp == "magna":
			run_command(['skewer', '-t', threads, '-q', qual_threshold, '-u', '-x', adapter_file,
				'-y', adapter_file_rev, p0_in, p1_in])

		# Rename trimmed files
		p1_new_path= os.path.join(dirs['trim'], p[0].replace('.fq', '.trimmed.fq'))
		p2_new_path= os.path.join(dirs['trim'], p[0].replace('1.sort.rePhred.fq',
			'2.sort.rePhred.trimmed.fq'))

		p1_trimmed = os.path.join(dirs['trim'], p[0].replace('.fq', '-trimmed-pair1.fastq'))
		p2_trimmed = os.path.join(dirs['trim'], p[0].replace('.fq', '-trimmed-pair2.fastq'))
		shutil.move(p1_trimmed, p1_new_path)
		shutil.move(p2_trimmed, p2_new_path)
		os.remove(p0_in)
		os.remove(p1_in)

	## Move Skewer logfiles
	for f in glob(os.path.join(dirs['trim'], '*.log')):
		in_log = os.path.join(dirs['trim'], os.path.basename(f))
		out_log = os.path.join(dirs['skewerlog'], os.path.basename(f))
		shutil.move(in_log, out_log)

	## Fastqc
	for f in os.listdir(dirs['trim']):
		if f.endswith('.fq'):
			in_fq = os.path.join(dirs['trim'], f)
			run_command(['fastqc', '-t', threads, '-o', dirs['trim_fastqc'], in_fq])

	run_comp_map()

def run_comp_map():
	clear_dirs([dirs['comp_map'], dirs['comp_maplog'], dirs['complog'],
			dirs['mtdna_complog'], dirs['mtdna_comp_map'], dirs['nuc_comp_map']])

	align(dirs['trim'], dirs['comp_map'], merged_ref_index, dirs['comp_maplog'])
	to_bam_and_sort(dirs['comp_map'])

	bams = glob(dirs['comp_map']+'/*.bam')

	run_map_to_og_and_rot(True)

def run_map_to_og_and_rot(skip_extract):
	clear_dirs([dirs['map_rot_og'], dirs['og_map'], dirs['oglog'],
        dirs['og_maplog'], dirs['og_mrkduplog'], dirs['rot_map'],
        dirs['rot_maplog'], dirs['rot_mrkduplog'],
		dirs['rotlog']])

	if spp == 'pulex':
		mit = '\"gi|5835848|ref|NC_000844.1|\"'
	elif spp == 'magna':
		mit = '\"NC_026914.1\"'

	if skip_extract:

		## Extract reads mapped to mitochondria
		for f in os.listdir(dirs['comp_map']):
			if f.endswith('.bam'):
				in_bam = os.path.join(dirs['comp_map'], f)
				out_bam = os.path.join(dirs['mtdna_comp_map'],
						f.replace('.bam', '.mapped_to_mtdna.bam'))
				os.system('samtools view -b -@ '+threads+' '+in_bam+' '+mit+' > '+out_bam)

				run_command(['samtools', 'index', out_bam])
		## Extract reads mapped to nuclear genome
		for f in os.listdir(dirs['comp_map']):
			if f.endswith('.bam'):
				in_bam = os.path.join(dirs['comp_map'], f)
				out_bam = os.path.join(dirs['nuc_comp_map'], f.replace('.bam', '.mapped_to_nuc.bam'))
				os.system("samtools idxstats "+in_bam+" | cut -f 1 | grep -v "+mit+ \
						" | xargs samtools view -b -@ "+threads+" "+in_bam+" > "+out_bam)
				run_command(['samtools', 'index', out_bam])


		for f in os.listdir(dirs['mtdna_comp_map']):
			if f.endswith('.bam'):
				print(f)
				bam_to_fastq(os.path.join(dirs['mtdna_comp_map'], f), dirs['mtdna_comp_map'])

	align(dirs['mtdna_comp_map'], dirs['og_map'], mtdna_og_index, dirs['og_maplog'])
	align(dirs['mtdna_comp_map'], dirs['rot_map'], mtdna_rot_index, dirs['rot_maplog'])
	to_bam_and_sort(dirs['og_map'])
	post_map(dirs['og_map'], dirs['og_mrkduplog'], mtdna_og_ref)
	to_bam_and_sort(dirs['rot_map'])
	post_map(dirs['rot_map'], dirs['rot_mrkduplog'], mtdna_rot_ref)

	run_call_variants()

def run_call_variants():
	clear_dirs([dirs['og_var'], dirs['rot_var'], dirs['og_var_unique'],
		dirs['rot_var_unique']])

	call_variants(dirs['og_map'], mtdna_og_ref, dirs['og_var'])
	call_variants(dirs['rot_map'], mtdna_rot_ref, dirs['rot_var'])


if __name__ == "__main__":
    arg = docopt(__doc__)


	# Get species parameter
	if arg['<species>' == "pulex" or arg['<species>' == "magna":
		spp = arg['<species>']
	else:
		print('Species should be one of pulex or magna')
		sys.exit()

	# Get isolate parameter
	if arg['isolate'] in ['LIN', 'TCO', 'I', 'G', 'F']:
		isolate = arg['isolate']
	else:
		print('Isolate should be one of (LIN/TCO/I/G/F).')
		sys.exit()

	print("Running with the following parameters-")
	print("Species: "+spp)
	print("Isolate: "+isolate)

	## Basic directory locations
	dirs = {}
	# Home directory for the whole MA project
	homedir = os.path.join(os.path.expanduser('~'),
                        'fmacrae/daphnia_mitochondria_MA/data')
	print('Project root: '+homedir)

	# Directory containing unprocessed read files for isolate
	if spp == 'pulex':
		if isolate == 'LIN':
			#dirs['root'] = os.path.join(homedir, 'illumina_pulex_indiana/LIN_reads')
			dirs['root'] = os.path.join(homedir, 'illumina_pulex_indiana/LIN_reads')
		elif isolate == 'TCO':
			dirs['root'] = os.path.join(homedir, 'illumina_pulex_indiana/TCO_reads')
	elif spp == 'magna':
		if isolate == 'I':
			dirs['root'] = os.path.join(homedir, 'illumina_magna_150bp/I')
		elif isolate == 'F':
			dirs['root'] = os.path.join(homedir, 'illumina_magna_150bp/F')
		elif isolate == 'G':
			dirs['root'] = os.path.join(homedir, 'illumina_magna_150bp/G')

	# Reference Sequence Locations
	if spp == "pulex":
		mtdna_app_ref = os.path.join(homedir, 'ref/pulex/mit/pulex_mtdna.app.fa')
		mtdna_app_index = os.path.join(homedir, 'ref/pulex/mit/pulex_mtdna.app')
		mtdna_og_ref = os.path.join(homedir, 'ref/pulex/mit/pulex_mtdna.og.fa')
		mtdna_og_index = os.path.join(homedir, 'ref/pulex/mit/pulex_mtdna.og')
		mtdna_rot_ref = os.path.join(homedir, 'ref/pulex/mit/pulex_mtdna.rot.fa')
		mtdna_rot_index = os.path.join(homedir, 'ref/pulex/mit/pulex_mtdna.rot')
		nuc_ref = os.path.join(homedir, 'ref/pulex/nuc/pulex_scaffolds.fa')
		nuc_ref_index = os.path.join(homedir, 'ref/pulex/nuc/pulex_scaffolds')
		merged_ref = os.path.join(homedir, 'ref/pulex/nuc_mit_merged/nuc_and_mito.fa')
		merged_ref_index = os.path.join(homedir, 'ref/pulex/nuc_mit_merged/nuc_and_mito')
	elif spp == "magna":
		mtdna_app_ref = os.path.join(homedir, 'ref/magna/mit/magna_mtdna.app.fa')
		mtdna_app_index = os.path.join(homedir, 'ref/magna/mit/magna_mtdna.app')
		mtdna_og_ref = os.path.join(homedir, 'ref/magna/mit/magna_mtdna.og.fa')
		mtdna_og_index = os.path.join(homedir, 'ref/magna/mit/magna_mtdna.og')
		mtdna_rot_ref = os.path.join(homedir, 'ref/magna/mit/magna_mtdna.rot.fa')
		mtdna_rot_index = os.path.join(homedir, 'ref/magna/mit/magna_mtdna.rot')
		nuc_ref = os.path.join(homedir, 'ref/magna/nuc/dmagna-v2.4-20100422-assembly.fa')
		nuc_ref_index = os.path.join(homedir, 'ref/magna/nuc/dmagna-v2.4-20100422-assembly')
		merged_ref = os.path.join(homedir, 'ref/magna/nuc_mit_merged/nuc_and_mito.fa')
		merged_ref_index = os.path.join(homedir, 'ref/magna/nuc_mit_merged/nuc_and_mito')

	# Adapter sequence locations - Only specified for D. magna. Default settings work for pulex.
	if spp == "magna":
		adapter_file = os.path.join(homedir, 'adapters/magna/adapters_short.fa')
		adapter_file_rev = os.path.join(homedir, 'adapters/magna/adapters_rc_short.fa')

	# Processing tree directory locations
	# QC
	dirs['trim'] = os.path.join(dirs['root'], 'processed')
	dirs['trimlog'] = os.path.join(dirs['trim'], 'logs')
	dirs['skewerlog'] = os.path.join(dirs['trimlog'], 'skewer')
	dirs['trim_fastqc'] = os.path.join(dirs['trim'], 'fastqc')

	# Competitive alignment and variant calling
	dirs['comp_map'] = os.path.join(dirs['trim'], 'competitive_map')
	dirs['complog'] = os.path.join(dirs['comp_map'], 'logs')
	dirs['comp_maplog'] = os.path.join(dirs['complog'], 'map')
	dirs['mtdna_comp_map'] = os.path.join(dirs['comp_map'], 'mtdna')
	dirs['mtdna_complog'] = os.path.join(dirs['mtdna_comp_map'], 'logs')
	dirs['map_rot_og'] = os.path.join(dirs['mtdna_comp_map'], 'map_to_rot_and_og_ref')
	dirs['og_map'] = os.path.join(dirs['map_rot_og'], 'og')
	dirs['oglog'] = os.path.join(dirs['og_map'], 'logs')
	dirs['og_maplog'] = os.path.join(dirs['oglog'], 'map')
	dirs['og_mrkduplog'] = os.path.join(dirs['oglog'], 'mrkdup')
	dirs['rot_map'] = os.path.join(dirs['map_rot_og'], 'rot')
	dirs['rotlog'] = os.path.join(dirs['rot_map'], 'logs')
	dirs['rot_maplog'] = os.path.join(dirs['rotlog'], 'map')
	dirs['rot_mrkduplog'] = os.path.join(dirs['rotlog'], 'mrkdup')
	dirs['nuc_comp_map'] = os.path.join(dirs['comp_map'], 'nuc')
	dirs['og_var'] = os.path.join(dirs['og_map'], 'variants')
	dirs['rot_var'] = os.path.join(dirs['rot_map'], 'variants')
	dirs['og_var_unique'] = os.path.join(dirs['og_var'], 'unique')
	dirs['rot_var_unique'] = os.path.join(dirs['rot_var'], 'unique')

	## Sequence analysis parameters

	# Num threads = num cores
	threads = str(multiprocessing.cpu_count())

	# Threshold for quality trimming
	qual_threshold = '30'

	# Bowtie mapping preset
	bowtie2_preset = '--very-sensitive'

	#Whether to remove duplicates or just mark them
	removeDups = 'True'

	# Region to call variants in og and rotated references
	region = os.path.join(homedir, "bed/"+spp+"_call_region.bed")
	init_dirs(dirs)
	startdir = os.getcwd()

	run_pipeline(startdir)


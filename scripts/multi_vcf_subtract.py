#!/usr/bin/env python

## Given a list of VCF.gz files, repeatedly call lofreq subtract to find the 
## variants which are unique to each sample.
import sys
from os import getcwd, path, system, makedirs, remove
from shutil import copy2, move

# Return the sample name of a file
def sample_name(f):
	return path.basename(f).split('.')[0]

def extract_vcf_header(f):
	header = ""
	with open(f, 'r') as vcf_file:
		for line in vcf_file:
			if line.startswith("#"):
				header = header + line
	
	return header

def prepend_to_file(string, filename):
	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0,0)
		f.write(string.rstrip('\r\n') + '\n' + content)

if __name__ == "__main__":

	cwd = getcwd()

	inFiles = []

	print("Finding unique variants in:")
	for i in range(1, len(sys.argv)):
		basePath = path.basename(sys.argv[i])
		print(basePath)
		fullPath = path.abspath(sys.argv[i])
		inFiles.append(fullPath)

	n = len(inFiles)

	odir = path.join(cwd, "unique")
	if not path.exists(odir):
		makedirs(odir)

	inZipped = []
	for f in inFiles:
		bgzipOut = f+".gz"

		if path.exists(bgzipOut):
			remove(bgzipOut)
			remove(bgzipOut+".tbi")
		system("bgzip -c "+f+" > "+bgzipOut)
		system("tabix "+bgzipOut)

		inZipped.append(bgzipOut)

	inZipStr = " ".join(inZipped)
	sharedSites = path.join(odir, "shared.vcf")
	system("bcftools isec --nfiles="+str(n)+" -o "+sharedSites+" -w 2 "+inZipStr)

	for f in inFiles:
		header = extract_vcf_header(f)
		unique = path.join(odir, sample_name(f)+".unique.vcf")
		system("bedtools subtract -a "+f+" -b "+sharedSites+" > "+unique)
		prepend_to_file(header, unique)
	

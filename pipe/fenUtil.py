import glob, subprocess, os, shutil

## Run a bash command as a subprocess
def run_command(command):
	ps = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = ps.communicate()
	print(out)
	print(err)
	
	return out, err

## Returns a list of tuples of paired files in the current directory 
## Catered to daphnia pulex and megadaph data (doesn't work with arbitrary 
## naming conventions)
def pair_files():
	## Grab an arbitrary file to assert filename structure
	all_files = glob.glob('*.fq*') + glob.glob("*.fastq*")
	first = all_files[0]
	
	## Pulex data naming convention
	if "_1." in first or "_2." in first:
        	p1 = sorted(glob.glob('*_1.*'))
        	p2 = sorted(glob.glob('*_2.*'))
        	return zip(p1, p2)
	
	## Megadaph naming convention
	elif "lane" in first:
        	p1 = sorted(glob.glob('*_R1_*'))
        	p2 = sorted(glob.glob('*_R2_*'))
        	return zip(p1, p2)
		
## Delete output from a previous analysis and initialize a new empty directory
def init_dirs(outdir):
        if os.path.exists(outdir):
                shutil.rmtree(outdir)

        os.makedirs(outdir)

## Return the library name of a given file
def library_name(f):
	return f.split('_', 1)[0]

## Return the sample name of a given file
def sample_name(f):
	return f.split('.', 1)[0][-1:]


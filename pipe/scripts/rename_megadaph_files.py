#!/usr/bin/env python
import os

if __name__ == "__main__":
	for f in os.listdir(os.getcwd()):
		if f.endswith(".fq") or f.endswith(".fastq") or f.endswith(".html"):
			in_fq = os.path.join(os.getcwd(), f)
			
			hyph_split = f.split('-')[-1]
			if "_EC" in f or "_SC" in f:
				underscore_split = hyph_split.split('_')[2:4]
				sample = underscore_split[0]+"_"+underscore_split[1]
			else:
				sample = hyph_split.split('_')[2]

			if "_R1_" in f:
				orient = "1"
			elif "_R2_" in f:
				orient = "2"
			
			if f.endswith(".fq") or f.endswith(".fastq"):
				out_fq = os.path.join(os.getcwd(), sample+"_"+orient+".fq")
			else:
				out_fq = os.path.join(os.getcwd(), sample+"_"+orient+".html")
			os.rename(in_fq, out_fq)



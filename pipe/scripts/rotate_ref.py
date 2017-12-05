#!/usr/bin/env python

## Assumes sequence is not split across multiple lines

from __future__ import print_function
import sys

def strip_new_lines(lst):
    return map(lambda s: s.strip(), lst)

if __name__ == "__main__":

    fname = sys.argv[1]

    # Read input
    with open(fname, 'r') as f:
        content = f.readlines()
        content = strip_new_lines(content)

    header = content[0]
    seq = ''.join(content[1:])
    bp = len(seq)
    midpoint = int(bp/2)

    # Rotate sequence
    first_half = seq[0:midpoint]
    second_half = seq[midpoint:]
    rot_seq = second_half + first_half
    header = header + " " + str(midpoint) + ":" + str(bp) + " rot to start"

    # Rename out file
    if fname.endswith(".fasta"):
        rotname = fname[:-5] + "rot.fasta"
    elif fname.endswith(".fa"):
        rotname = fname[:-2] + "rot.fa"

    # Write out file
    with open(rotname, 'w') as f:
        print(header, file=f)
        print(rot_seq, file=f)

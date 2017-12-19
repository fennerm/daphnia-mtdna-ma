#!/usr/bin/env python
## Assumes sequence is not split across multiple lines

from __future__ import print_function
import sys

from plumbum import (
    FG,
    local,
    )

SCRIPTDIR = local.path("scripts")

def strip_new_lines(lst):
    return [s.strip() for s in lst]

def main(fname, rotname):

    # Read input
    with open(fname, 'r') as f:
        content = f.readlines()
        content = strip_new_lines(content)

    header = content[0][1:]

    seq = ''.join(content[1:])
    bp = len(seq)
    midpoint = int(bp/2)

    # Rotate sequence
    first_half = seq[0:midpoint]
    second_half = seq[midpoint:]
    rot_seq = second_half + first_half
    old_header = header.partition(' ')[0]
    header = ''.join([
        '>', str(midpoint), ":", str(bp), "_rot_to_start", '_',
        old_header])

    # Write out file
    with open(rotname, 'w') as f:
        print(header, file=f)
        print(rot_seq, file=f)


if __name__ == "__main__":
    fname = sys.argv[1]
    rotname = sys.argv[2]
    main(fname, rotname)


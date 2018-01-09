#!/usr/bin/env python
## Assumes sequence is not split across multiple lines

from __future__ import print_function
from math import ceil
import sys

from plumbum import (
    FG,
    local,
    )

SCRIPTDIR = local.path("scripts")

def strip_new_lines(lst):
    return [s.strip() for s in lst]

def main(fname, rotname = None):
    # If no rotname given, output is written to stdout

    # Read input
    with open(fname, 'r') as f:
        content = f.readlines()
        content = strip_new_lines(content)

    header = content[0][1:]

    seq = ''.join(content[1:])
    bp = len(seq)
    midpoint = int(ceil(bp/2))

    # Rotate sequence
    first_half = seq[0:midpoint]
    second_half = seq[midpoint:]
    rot_seq = second_half + first_half
    old_header = header.partition(' ')[0]

    # We convert to 1 based indexing
    header = ''.join([
        '>', str(midpoint + 1), ":", str(bp), "_rot_to_start", '_',
        old_header])

    # Write out file

    if rotname:
        with open(rotname, 'w') as f:
            print(header, file=f)
            print(rot_seq, file=f)
    else:
        sys.stdout.write(header + '\n')
        sys.stdout.write(rot_seq + '\n')


if __name__ == "__main__":
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        main(sys.argv[1])

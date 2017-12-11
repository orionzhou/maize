#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'clean sequence IDs in a fasta file'
    )
    parser.add_argument(
        'fi', help = 'input file (.fas)'
    )
    parser.add_argument(
        'fo', help = 'output file (.fas)'
    )
    args = parser.parse_args()
    
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith(">"):
            fho.write(line.rstrip(":.\n")+"\n")
        else:
            fho.write(line)
    fhi.close()
    fho.close()

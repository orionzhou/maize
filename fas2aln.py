#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import AlignIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert alignment file in fasta format to clustal format'
    )
    parser.add_argument(
        'fi', help = 'input alignment (.fas)'
    )
    parser.add_argument(
        'fo', help = 'output alignment (.aln)'
    )
    args = parser.parse_args()

    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")

    alns = AlignIO.parse(fhi, "fasta")
    AlignIO.write(alns, fho, "clustal")
    fhi.close()
    fho.close()

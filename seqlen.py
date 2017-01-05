#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'report fasta sequence lengths'
    )
    parser.add_argument(
        'fi', nargs = '?', help = 'input file (fasta)'
    )
    parser.add_argument(
        'fo', nargs = '?', help = 'output file (tbl)'
    )
    args = parser.parse_args()

    if args.fi and args.fi != '-' and args.fi != 'stdin':
        sys.stdin = open(args.fi, "r")
    if args.fo and args.fo != '-' and args.fo != 'stdout':
        sys.stdout = open(args.fo, "w")

    for seq in SeqIO.parse(sys.stdin, "fasta") :
        print("\t".join([seq.id, str(len(seq.seq))]))

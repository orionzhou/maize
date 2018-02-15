#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

def summary(args):
    if args.fi and args.fi != '-' and args.fi != 'stdin':
        sys.stdin = open(args.fi, "r")
    if args.fo and args.fo != '-' and args.fo != 'stdout':
        sys.stdout = open(args.fo, "w")
    
    print("seqid\tsize\tdesc")
    for seq in SeqIO.parse(sys.stdin, "fasta") :
        print("\t".join([seq.id, str(len(seq.seq)), seq.description]))

def fas2aln(args):
    from Bio import AlignIO
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    alns = AlignIO.parse(fhi, "fasta")
    AlignIO.write(alns, fho, "clustal")
    fhi.close()
    fho.close()

def rmdot(args):
    from string import maketrans
    tt = maketrans(".", "-")
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith('>'):
            fho.write(line)
        else:
            fho.write(line.translate(tt))
    fhi.close()
    fho.close()

def cleanid(args):
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith(">"):
            fho.write(line.rstrip(":.\n")+"\n")
        else:
            fho.write(line)
    fhi.close()
    fho.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            prog = "python -m maize.formats.fasta",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fasta utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("summary", help = "Report length and description for each sequence")
    sp1.add_argument('fi', nargs = '?', help = 'input file (fasta)')
    sp1.add_argument('fo', nargs = '?', help = 'output file (tsv)')
    sp1.set_defaults(func = summary)

    sp11 = sp.add_parser("fas2aln", help = 'convert fasta alignment file to clustal format')
    sp11.add_argument('fi', help = 'input alignment (.fas)')
    sp11.add_argument('fo', help = 'output alignment (.aln)')
    sp11.set_defaults(func = fas2aln)

    sp31 = sp.add_parser("rmdot", help = 'replace periods (.) in an alignment fasta by dashes (-)')
    sp31.add_argument('fi', help = 'input fasta file')
    sp31.add_argument('fo', help = 'output fasta file')
    sp31.set_defaults(func = rmdot)
    
    sp32 = sp.add_parser("cleanid", help = 'clean sequence IDs in a fasta file')
    sp32.add_argument('fi', help = 'input fasta file')
    sp32.add_argument('fo', help = 'output fasta file')
    sp32.set_defaults(func = cleanid)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


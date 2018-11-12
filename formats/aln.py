#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import AlignIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert clustal (.aln) file to other formats'
    )
    parser.add_argument(
        'fi', help = 'input alignment (.aln)'
    )
    parser.add_argument(
        'fo', help = 'output alignment file'
    )
    parser.add_argument(
        '--len', dest='len', type=int, default=10, help='id length (def: 10)'
    )
    parser.add_argument(
        '--fmt', dest='fmt', type=int, default=1, help='1: fasta (defalt); 2: phylip; 3: nexus; 4: stockholm'
    )
    args = parser.parse_args()

    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    fmt = args.fmt

    alns = AlignIO.parse(fhi, "clustal")
    if fmt == 1: 
       AlignIO.write(alns, fho, "fasta")
    elif fmt == 2:
       AlignIO.write(alns, fho, "phylip-relaxed")
    elif fmt == 3:
       AlignIO.write(alns, fho, "nexus")
    elif fmt == 4:
       AlignIO.write(alns, fho, "stockholm")
    else:
        print("unknown fmt: %s" % fmt)
        sys.exit(1)

    fhi.close()
    fho.close()

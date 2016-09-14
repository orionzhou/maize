#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import Phylo

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert alignment file in fasta format to clustal format'
    )
    parser.add_argument(
        'fi', help = 'input phylogeny file(.nwk)'
    )
    parser.add_argument(
        'fo', help = 'output phylogeny (.nex)'
    )
    args = parser.parse_args()

    tree = Phylo.read(args.fi, "newick")
    Phylo.write(tree, args.fo, "nexus")

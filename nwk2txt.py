#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import Phylo

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'extract leaf(tip) labels from newick tree file'
    )
    parser.add_argument(
        'fi', nargs = '?', help = 'input file (newick)'
    )
    parser.add_argument(
        'fo', nargs = '?', help = 'output file (txt)'
    )
    args = parser.parse_args()

    (fi, fo) = (args.fi, args.fo)
    tree = Phylo.read(fi, "newick")
    
    labels = []
    for leaf in tree.get_terminals():
        labels.append(leaf.name)
    fho = open(fo, "w")
    for label in reversed(labels):
        print >>fho, label
    fho.close()

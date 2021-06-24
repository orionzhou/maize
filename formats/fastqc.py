#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read and write FASTQC output.

Library dependency: Fadapa
"""

import os.path as op
import sys
import logging

from jcvi.apps.base import sh, mkdir
from fadapa import Fadapa

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read and write fastqc output'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('size_dist', help='write read length distribution',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'fastqc data file ("fastqc_data.txt")')
    sp1.add_argument('fo', help = 'output tsv file')
    sp1.set_defaults(func = size_dist)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

def size_dist(args):
    fi = Fadapa(args.fi)
    fho = open(args.fo, "w")
    for r in fi.clean_data("Sequence Length Distribution"):
        fho.write("\t".join(r) + "\n")
    fho.close()

if __name__ == '__main__':
    main()

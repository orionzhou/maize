#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def lookup(args):
    genome = args.genome
    if args.genome.upper() == ['B73','W22','PH207']:
        genome = "Zmays_" + genome
    elif args.genome.upper() == 'MO17':
        genome = 'Zmays_Mo17'
    fi = f"{args.dirg}/data2/loci/{genome}.tsv"
    sh(f"grep {args.word} {fi}", log=False)


if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "look up locus using GeneID, name, or description"
    )

    ps.add_argument('word', help = 'work to look up')
    ps.add_argument('--genome', '-g', default='Zmays_B73', help = 'genome')
    ps.add_argument('--fmt', default='long', help = 'output format')
    ps.add_argument('--dirg', default='/home/springer/zhoux379/projects/genome', help = 'genome directory')

    args = ps.parse_args()
    lookup(args)


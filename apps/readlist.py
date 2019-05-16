#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Readlist utilities
"""

import os.path as op
import sys
import re
import logging
import pandas as pd

from jcvi.formats.fastq import iter_fastq
from jcvi.utils.cbook import SummaryStats
from maize.formats.base import must_open, DictFile

def add_stat(args):
    cvt = {k: int for k in 'Replicate'.split()}
    sl = pd.read_csv(args.fi, sep="\t", header=0, converters=cvt)
    firstN = 10000

    sl['spots'] = [0] * len(sl.index)
    sl['avgLength'] = [0] * len(sl.index)
    for i in range(len(sl)):
        sid = sl['SampleID'][i]
        fq = ''
        if sl['paired'][i]:
            r1, r2 = sl['r1'][i], sl['r2'][i]
            fq = r1
        else:
            fq = sl['r0'][i]

        nrcd = 0
        L = []
        for rec in iter_fastq(fq):
            if not rec:
                break
            nrcd += 1
            if nrcd <= firstN:
                L.append(len(rec))

        avgLength = SummaryStats(L).mean
        if sl['paired'][i]:
            avgLength = avgLength * 2

        print("\t".join(str(x) for x in (sid, nrcd, avgLength)))
        sl.at[i, 'spots'] = nrcd
        sl.at[i, 'avgLength'] = avgLength

    sl.to_csv(args.fo, sep="\t", header=True, index=False)

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fastq utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('add_stat', help='add fastq stats [spots, avgLength]',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input read list')
    sp1.add_argument('fo', help = 'output read list with stats added')
    sp1.set_defaults(func = add_stat)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

if __name__ == '__main__':
    main()


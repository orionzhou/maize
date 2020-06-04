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

from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles

def venn3_coord(args):
    fhi = open(args.fi, 'r')
    s1 = fhi.readline().strip().split(",")
    s2 = fhi.readline().strip().split(",")
    s3 = fhi.readline().strip().split(",")
    fhi.close()
    s1, s2, s3 = set(s1), set(s2), set(s3)
    v = venn3([s1, s2, s3], ('A','B','C'))
    fho1 = open(args.fo1, 'w')
    for xy, l in zip(v.centers, v.radii):
        x, y = xy
        fho1.write("%s\t%s\t%s\n" % (x, y, l))
    fho1.close()
    fho2 = open(args.fo2, 'w')
    for xyl in v.subset_labels:
        x, y = xyl.get_position()
        l = xyl.get_text()
        fho2.write("%s\t%s\t%s\n" % (x, y, l))
    fho2.close()

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
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = '3-way venn-diagram'
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('coord', help='compute venn3 coordinates',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input file containing sets')
    sp1.add_argument('fo1', help = 'output circle coordinates')
    sp1.add_argument('fo2', help = 'output label coordinates')
    sp1.set_defaults(func = venn3_coord)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

if __name__ == '__main__':
    main()


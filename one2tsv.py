#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = 'format conversion'
    )
    parser.add_argument(
            'fi', help = 'input file (cluster-one output in *.csv)'
    )
    parser.add_argument(
            'fo', help = 'output file (*.tsv)'
    )
    parser.add_argument(
            '--maxp', type=float, default=0.05, help='P-value threshold (default: 0.1)'
    )
    args = parser.parse_args()
    
    if args.fi and args.fi != '-' and args.fi != 'stdin':
        sys.stdin = open(args.fi, "r")
    if args.fo and args.fo != '-' and args.fo != 'stdout':
        sys.stdout = open(args.fo, "w")
     
    print("grp\tgid")
    grp = 1
    for line in sys.stdin:
        ps = line.strip("\n").split(",")
        if ps[0] == 'Cluster':
            continue
        mid, size, density, iwt, ewt, quality, pval, gidstr = ps
        gids = gidstr.replace("\"", "").split(" ")
        if float(pval) >= args.maxp or len(gids) < 5:
            continue
        for gid in gids:
            print("%d\t%s" % (grp, gid))
        grp += 1

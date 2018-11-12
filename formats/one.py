#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cluter-One format
"""

import os
import os.path as op
import sys
import logging
import re

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

def one2tsv(args):
    fhi = must_open(args.mcl)
    print("grp\tgid")
    grp = 1
    for line in fhi:
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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Cluster-One format utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv", help = "one -> tsv")
    sp1.add_argument('fi', help = 'input file (cluster-one output in *.csv)')
    sp1.add_argument('--maxp', type=float, default=.05, help='p-value thresh')
    sp1.set_defaults(func = one2tsv)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



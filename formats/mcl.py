#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MCL format
"""

import os
import os.path as op
import sys
import logging
import re

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

def mcl2tsv(args):
    fhi = must_open(args.mcl)
    print("grp\tgid")
    grp = 1
    for line in fhi:
        line = line.strip("\n")
        gids = line.split("\t")
        if len(gids) < 5:
            continue
        for gid in gids:
            print("%d\t%s" % (grp, gid))
        grp += 1

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'MCL format utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv", help = "mcl -> tsv")
    sp1.add_argument('mcl', help = 'input MCL file')
    sp1.set_defaults(func = mcl2tsv)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



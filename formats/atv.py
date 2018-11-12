#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ATV (alignment-tsv) format
"""

import os
import os.path as op
import sys
import logging
import re

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

def filter(args):
    fhi = must_open(args.fi)
    line = fhi.readline()
    print(line.strip("\n"))
    pqid = ''
    pscore = 0
    lines = []
    for line in sys.stdin:
        line = line.strip("\n")
        qName, qStart, qEnd, qSize, strand, \
        tName, tStart, tEnd, tSize, \
        alnLen, match, misMatch, baseN, \
        qNumIns, tNumIns, qBaseIns, tBaseIns, \
        ident, score, qLoc, tLoc = line.split("\t")
        #print(line)
        #print(qSize)
        if float(ident) < args.ident: continue
        if int(alnLen)/int(qSize) < args.cov: continue
        if int(match) < args.match: continue
        if pqid == '':
            pqid = qName
            pscore = score
            lines.append(line)
        elif qName != pqid:
            print("\n".join(lines))
            pqid = qName
            pscore = score
            lines = [line]
        else:
            if args.best:
                if score > pscore:
                    lines = [line]
                    pscore = score
                elif score == pscore:
                    lines.append(line)
            else:
                lines.append(line)
    print("\n".join(lines))
 
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'ATV (alignment-tsv) format utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("filter", help = "filter an alignment file",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input alignment file (tsv)')
    sp1.add_argument('--match', type=int, default=1, help='min matches')
    sp1.add_argument('--ident', type=float, default=0.8, help='min identity')
    sp1.add_argument('--cov', type=float, default=0.8, help='min query coverage')
    sp1.add_argument('--best', action = 'store_true', help='only keep best hit(s)')
    sp1.set_defaults(func = filter)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

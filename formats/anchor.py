#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import argparse
from maize.formats.base import LineFile, must_open, is_number, get_number, ndigit, prettysize
from jcvi.compara.synteny import AnchorFile

def anchor2tsv(args):
    anchors = AnchorFile(args.fi)
    blocks = anchors.blocks
    i = 1
    fmt = "b%%0%dd" % ndigit(len(blocks))
    print("\t".join('bid gid1 gid2 score'.split()))
    for block in blocks:
        bid = fmt % i
        for line in block:
            a, b, score = line
            print("\t".join((bid, a, b, score)))
        i += 1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'anchor file utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv", help = "convert anchor file to *.tsv")
    sp1.add_argument('fi', help = 'input *.anchor file')
    sp1.set_defaults(func = anchor2tsv)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



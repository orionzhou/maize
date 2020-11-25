#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re

from jcvi.formats.base import must_open, read_until, read_block

def dreme2tsv(args):
    fh = must_open(args.fi)

    start = fh.tell()
    pos_c, neg_c = '', ''
    for line in fh:
        ps = re.split(r" +", line)
        if len(ps) >= 2 and ps[1].startswith('positives'):
            pos_c = ps[2]
        elif len(ps) >= 2 and ps[1].startswith('negatives'):
            neg_c = ps[2]
            break
    fh.seek(start)

    i = 0
    fho = must_open(args.fo, 'w')
    fho.write("mid\tre\tseq\tseq_rc\tpos\tneg\tpos_c\tneg_c\tpval\teval\n")
    for head, lines in read_block(fh, 'MOTIF'):
        if not head: break
        i += 1
        mtf, conseq, mid = head.split(" ")
        for line in lines:
            if line.startswith("#"):
                line = line.replace("#", '').strip()
                ps = re.split(r" +", line)
                if ps[0] in ['Word', 'Stopping', 'Running']: continue
                tag = ''
                if ps[0] == 'BEST':
                    tag = 'RE'
                    ps = ps[1:]
                seq, seq_rc, pos, neg, pval, eval = ps
                fho.write("\t".join([mid, tag, seq, seq_rc, pos, neg, pos_c, neg_c, pval, eval])+'\n')
    fh.close()

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'dreme utilities'
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv", help = "dreme -> tsv")
    sp1.add_argument('fi', help = 'input dreme file')
    sp1.add_argument('fo', help = 'output tsv file')
    sp1.set_defaults(func = dreme2tsv)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()



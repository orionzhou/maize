#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open
from maize.utils.location import locAry2Str, locStr2Ary

def gtb2tsv(args):
    fhi = must_open(args.fi)
    print("\t".join("gid tid ttype etype chrom start end srd fam note".split()))
    for line in fhi:
        line = line.strip("\n")
        if line.startswith("#") or line.startswith("id"):
            continue
        ary = line.split("\t")
        if len(ary) < 18:
            print("less than 18 columns:\n%s" % line)
            continue
        tid, gid, seqid, tbeg, tend, srd, \
                locES, locIS, locCS, loc5S, loc3S, phase, \
                src, conf, cat1, cat2, cat3, note = ary
        tbeg, tend = int(tbeg), int(tend)
        if cat1 == 'mRNA':
            assert locCS, "no CDS for %d" % tid
        else:
            assert locES, "no exon for %d" % tid
        ldic = { 'exon': locES, 'cds': locCS, \
                'utr5': loc5S, 'utr3': loc3S, 'intron':locIS }
        for etype, locS in ldic.items():
            if not locS:
                continue
            for rbeg, rend in locStr2Ary(locS):
                beg, end = 0, 0
                if srd == "-":
                    beg, end = tend - rend + 1, tend - rbeg + 1
                else:
                    assert srd == '+', "unknown strand: %s for %s" % (srd, tid)
                    beg, end = tbeg + rbeg - 1, tbeg + rend - 1
                fields = [gid, tid, cat1, etype, seqid, str(beg), str(end), srd, cat3, note]
                print("\t".join(fields)) 
    fhi.close()

if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Gtb utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'gtb to tsv'
    )
    sp1.add_argument('fi', help = 'input Gtb file')
    sp1.set_defaults(func = gtb2tsv)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


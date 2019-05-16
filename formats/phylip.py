#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
from Bio import AlignIO

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open, ndigit, prettysize

def clean(args):
    reg = re.compile("[^ATCGN]")
    fh = must_open(args.fi)
    alns = AlignIO.read(fh, "phylip-relaxed")
    for rcd in alns:
        rcd.seq = reg.subn("N", str(rcd.seq).upper())[0]
        #rcd = SeqRecord(Seq(newseq), id = rcd.id)
    AlignIO.write(alns, sys.stdout, "phylip-relaxed")
    fh.close()

def phy2fa(args):
    fh = must_open(args.fi)
    alns = AlignIO.read(fh, "phylip")
    AlignIO.write(alns, sys.stdout, "fasta")
    fh.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'phylip utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("clean", help = "Remove non-ACGT characters by N")
    sp1.add_argument('fi', help = 'input file (phylip)')
    sp1.set_defaults(func = clean)

    sp1 = sp.add_parser("2fa", help = "convert to fasta format")
    sp1.add_argument('fi', help = 'input file (phylip)')
    sp1.set_defaults(func = phy2fa)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


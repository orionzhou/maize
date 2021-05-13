#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
MAF format specification:
<http://genome.ucsc.edu/FAQ/FAQformat#format5>
"""

import sys

from bx import interval_index_file
from bx.align import maf

from maize.formats.base import BaseFile
from jcvi.formats.maf import Maf
from maize.apps.base import need_update
from maize.apps.lastz import blastz_score_to_ncbi_expectation, \
            blastz_score_to_ncbi_bits

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = ''
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('bed', help='convert MAF to BED format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = bed)

    sp1 = sp.add_parser('blast', help='convert MAF to BLAST tabular format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = blast)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

def bed(args):
    """
    %prog bed maffiles > out.bed

    Convert a folder of maf alignments to the bed features
    then useful to check coverage, etc.
    """
    p = OptionParser(bed.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    flist = args
    prefix = flist[0].split(".")[0]

    j = 0
    for f in flist:
        reader = Maf(f).reader
        for rec in reader:
            a, b = rec.components

            for a, tag in zip((a, b), "ab"):
                name = "{0}_{1:07d}{2}".format(prefix, j, tag)
                print("\t".join(str(x) for x in (a.src, a.forward_strand_start, a.forward_strand_end, name)))

            j += 1

def alignment_details(a, b):
    nmatch = 0
    nmismatch = 0
    ngaps = 0

    assert len(a) == len(b)
    l = len(a)

    for i in range(l):
        if a[i] == b[i]:
            nmatch += 1
        elif a[i] == "-" or b[i] == "-":
            ngaps += 1
        else:
            nmismatch += 1

    pctid = 100. * nmatch / l
    return pctid, nmismatch, ngaps

def maf_to_blast8(f):
    reader = Maf(f).reader
    for rec in reader:
        a, b = rec.components
        query = a.src
        subject = b.src
        qstart = a.forward_strand_start
        qstop = a.forward_strand_end
        sstart = b.forward_strand_start
        sstop = b.forward_strand_end
        score = rec.score

        evalue = blastz_score_to_ncbi_expectation(score)
        score = blastz_score_to_ncbi_bits(score)
        evalue, score = "{0:.2g}".format(evalue), "{0:.1f}".format(score)
        hitlen = len(a.text)

        pctid, nmismatch, ngaps = alignment_details(a.text, b.text)
        print("\t".join(str(x) for x in (query, subject, pctid, hitlen,
            nmismatch, ngaps, qstart, qstop, sstart, sstop,
            evalue, score)))

def blast(args):
    '''
    %prog blast maffiles > out.blast

    From a folder of .maf files, generate .blast file with tabular format.
    '''
    p = OptionParser(blast.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(p.print_help())

    flist = args

    for f in flist:
        maf_to_blast8(f)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import math
import logging
import numpy as np

from jcvi.formats.bed import Bed, BedLine, BedpeLine, BedEvaluate, BedSummary

from jcvi.formats.base import LineFile, must_open, is_number, get_number
from jcvi.formats.sizes import Sizes
from jcvi.utils.cbook import SummaryStats, thousands, percentage
from jcvi.utils.grouper import Grouper
from jcvi.utils.range import Range, range_union, range_chain, \
            range_distance, range_intersect
from jcvi.apps.base import sh, need_update, mkdir

def make_window(beg, end, winsize, winstep):
    import math
    size = end - beg + 1
    nwin = math.ceil((size - winsize) / winstep) + 1
    nwin = max(1, nwin)

    mergelast = False
    if nwin > 1:
        size_lastwin = size - (nwin-1) * winstep
        if float(size_lastwin) / winstep < 0.3:
            mergelast = True
            #eprint(size_lastwin)

    wins = []
    for i in range(0, nwin):
        wbeg = beg + winstep * i
        wend = beg + winstep * i + winsize - 1
        wend = min(end, wend)
        if i == nwin - 2 and mergelast:
            wend = end
        elif i == nwin - 1 and mergelast:
            continue
        wins.append([wbeg, wend])
    return wins

def filter(args):
    """
    %prog filter bedfile

    Filter the bedfile to retain records between certain size range.
    """
    minsize, maxsize = args.minsize, args.maxsize
    minaccn = args.minaccn
    minscore = args.minscore
    total = []
    keep = []
    for row in must_open(args.fi):
        try:
            b = BedLine(row)
        except IndexError:
            logging.error(row.strip())
            continue
        span = b.span
        total.append(span)
        if not minsize <= span <= maxsize:
            continue
        if minaccn and int(b.accn) < minaccn:
            continue
        if minscore and int(b.score) < minscore:
            continue
        keep.append(span)
        print(b)
    logging.debug("Stats: {0} features kept.".\
        format(percentage(len(keep), len(total))))
    logging.debug("Stats: {0} bases kept.".\
        format(percentage(sum(keep), sum(total))))

def binpacking(args):
    import binpacking
    fi, fo, diro = args.bed, args.chain, args.outdir
    n = args.N
    #digit = ndigit(n)
    #fmt = "part.%%0%dd.fna" % digit
    fmt = f"{args.pre}%d.bed"
    if not op.exists(diro):
        mkdir(diro)
    else:
        sh("rm -rf %s/*" % diro)

    sdic = dict()
    i = 1
    fho = must_open(fo, 'w')
    for b in Bed(fi):
        rid = "s%d" % i
        sdic[rid] = [b.seqid, b.start-1, b.end, b.span]
        nid = "%s-%d-%d" % (b.seqid, b.start, b.end)
        fho.write("\t".join(str(x) for x in [nid, 0, b.span, '+',
                                             b.seqid, b.start-1, b.end, rid]) + "\n")
        i += 1
    fho.close()

    zdic = {k: v[3] for k, v in sdic.items()}
    bins = binpacking.to_constant_bin_number(zdic, n)

    sizes = []
    for j in range(n):
        sizes.append(np.sum(list(bins[j].values())))
        fl = op.join(diro, fmt % (j+1))
        fhl = must_open(fl, 'w')
        for rid in bins[j].keys():
            b = sdic[rid][:-1]
            fhl.write("\t".join(str(x) for x in b) + "\n")
        fhl.close()
    print("size range: %s - %s" % (sizes[0], sizes[n-1]))

def size(args):
    # sizes = [b.end - b.start + 1 for b in Bed(args.fi)]
    # sizes = np.sort(sizes)
    # total = np.sum(sizes)
    # min_sizes = " ".join(str(x) for x in sizes[0:3])
    # max_sizes = " ".join(str(x) for x in sizes[-3:])
    # logging.debug("size range: %s - %s, total size: %d" % (min_sizes, max_sizes, total))
    cmd = "bioawk -t '{sum += $3 - $2} END {print sum/1000000000}'"
    sh(cmd + " " + args.fi)

def makewindow(args):
    fhi = must_open(args.fi)
    for line in fhi:
        ps = line.strip().split("\t")
        seqid, start, end = ps[0], int(ps[1]) + 1, int(ps[2])
        assert start <= end, "start={0} end={1}".format(start, end)
        size = end - start + 1
        wins = make_window(start, end, args.size, args.step)
        for wbeg, wend in wins:
            print("%s\t%d\t%d" % (seqid, wbeg-1, wend))

def bed_sum(beds, seqid=None, unique=True):
    if seqid:
        ranges = [(x.seqid, x.start, x.end) for x in beds \
                    if x.seqid == seqid]
    else:
        ranges = [(x.seqid, x.start, x.end) for x in beds]

    unique_sum = range_union(ranges)
    raw_sum = sum(x.span for x in beds)
    return unique_sum if unique else raw_sum

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'bed utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("size", help = "report total size of features")
    sp1.add_argument('fi', help = 'input *.bed file')
    sp1.set_defaults(func = size)

    sp1 = sp.add_parser('filter', help='filter bedfile to retain records between size range',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input *.bed file')
    sp1.add_argument('--minsize', '-min', default=0, type=int, help='minimum feature length')
    sp1.add_argument('--maxsize', '-max', default=1000000000, type=int, help = 'maximum feature length')
    sp1.add_argument('--minaccn', type=int, help='minimum value of accn, useful to filter based on coverage')
    sp1.add_argument('--minscore', type=int, help='minimum score')
    sp1.set_defaults(func = filter)

    sp1 = sp.add_parser("makewindow", help = "make sliding windows with given size and step",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input *.bed file')
    sp1.add_argument('--size', '-w', default = 100, type = int, help = 'window size')
    sp1.add_argument('--step', '-s', default = 100, type = int, help = 'window step')
    sp1.set_defaults(func = makewindow)

    sp1 = sp.add_parser("binpacking", help = 'use bin packing library to distribute records evenly',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('bed', help = 'input intervals to split (*.bed)')
    sp1.add_argument('chain', help = 'output chain file (*.chain)')
    sp1.add_argument('outdir', help = 'output directory')
    sp1.add_argument('--N', type = int, default = 10, help = 'number pieces to split')
    sp1.add_argument('--pre', default='p', help = 'prefix of output fasta')
    sp1.set_defaults(func = binpacking)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


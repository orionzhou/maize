#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import math
import logging
import numpy as np

from collections import defaultdict
from itertools import groupby

from maize.formats.base import LineFile, must_open, is_number, get_number
#from maize.formats.sizes import Sizes
#from maize.utils.iter import pairwise
#from maize.utils.cbook import SummaryStats, thousands, percentage
#from maize.utils.grouper import Grouper
#from maize.utils.range import Range, range_union, range_chain, \
#            range_distance, range_intersect
from maize.apps.base import sh, need_update


class BedLine(object):
    """
    the Bed format supports more columns. we only need
    the first 4, but keep the information in 'extra'.
    """
    __slots__ = ("seqid", "start", "end", "accn",
                 "extra", "score", "strand", "args", "nargs")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.nargs = nargs = len(args)
        self.seqid = args[0]
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        assert self.start <= self.end, \
                "start={0} end={1}".format(self.start, self.end)
        self.extra = self.accn = self.score = self.strand = None

        if nargs > 3:
            self.accn = args[3]
        if nargs > 4:
            self.score = args[4]
        if nargs > 5:
            self.strand = args[5]
        if nargs > 6:
            self.extra = args[6:]

        self.args = args

    def __str__(self):
        args = [self.seqid, self.start - 1, self.end]
        if self.accn is not None:
            args += [self.accn]
        if self.score is not None:
            args += [self.score]
        if self.strand is not None:
            args += [self.strand]
        if self.extra is not None:
            args += self.extra

        s = "\t".join(str(x) for x in args)
        return s

    __repr__ = __str__

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def span(self):
        return self.end - self.start + 1

    @property
    def range(self):
        strand = self.strand or '+'
        return (self.seqid, self.start, self.end, strand)

    @property
    def tag(self):
        return "{0}:{1}-{2}".format(self.seqid, self.start, self.end)

    def reverse_complement(self, sizes):
        size = sizes.get_size(self.seqid)

        start = size - self.end + 1
        end = size - self.start + 1
        self.start, self.end = start, end
        assert self.start <= self.end, \
                "start={0} end={1}".format(self.start, self.end)

        if self.strand:
            strand = {'+': '-', '-': '+'}[self.strand]

    def gffline(self, type='match', source='default'):
        score = "." if not self.score or \
                (self.score and not is_number(self.score)) \
                else self.score
        strand = "." if not self.strand else self.strand
        row = "\t".join((self.seqid, source, type,
            str(self.start), str(self.end), score,
            strand, '.', 'ID=' + self.accn))
        return row

class Bed(LineFile):
    def __init__(self, filename=None, key=None, juncs=False, include=None):
        super(Bed, self).__init__(filename)

        if not filename:
            return

        for line in must_open(filename):
            if line[0] == "#" or (juncs and line.startswith('track name')):
                continue
            b = BedLine(line)
            if include and b.accn not in include:
                continue
            self.append(b)

    def add(self, row):
        self.append(BedLine(row))

    def print_to_file(self, filename="stdout", sorted=False):
        if sorted:
            self.sort(key=self.key)

        fw = must_open(filename, "w")
        for b in self:
            if b.start < 1:
                logging.error("Start < 1. Reset start for `{0}`.".format(b.accn))
                b.start = 1
            print >> fw, b
        fw.close()

    def sum(self, seqid=None, unique=True):
        return bed_sum(self, seqid=seqid, unique=unique)

    @property
    def seqids(self):
        return natsorted(set(b.seqid for b in self))

    @property
    def accns(self):
        return natsorted(set(b.accn for b in self))

    @property
    def order(self):
        # get the gene order given a Bed object
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

    @property
    def order_in_chr(self):
        # get the gene order on a particular seqid
        res = {}
        self.sort(key=self.nullkey)
        for seqid, beds in groupby(self, key=lambda x: x.seqid):
            for i, f in enumerate(beds):
                res[f.accn] = (seqid, i, f)
        return res

    @property
    def bp_in_chr(self):
        # get the bp position on a particular seqid
        res = {}
        self.sort(key=self.nullkey)
        for seqid, beds in groupby(self, key=lambda x: x.seqid):
            for i, f in enumerate(beds):
                res[f.accn] = (seqid, (f.start + f.end) / 2, f)
        return res

    @property
    def simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]

    @property
    def links(self):
        r = []
        for s, sb in self.sub_beds():
            for a, b in pairwise(sb):
                r.append(((a.accn, a.strand), (b.accn, b.strand)))
        return r

    def extract(self, seqid, start, end):
        # get all features within certain range
        for b in self:
            if b.seqid != seqid:
                continue
            if b.start < start or b.end > end:
                continue
            yield b

    def sub_bed(self, seqid):
        # get all the beds on one chromosome
        for b in self:
            if b.seqid == seqid:
                yield b

    def sub_beds(self):

        self.sort(key=self.nullkey)
        # get all the beds on all chromosomes, emitting one at a time
        for bs, sb in groupby(self, key=lambda x: x.seqid):
            yield bs, list(sb)

    def get_breaks(self):
        # get chromosome break positions
        simple_bed = self.simple_bed
        for seqid, ranks in groupby(simple_bed, key=lambda x: x[0]):
            ranks = list(ranks)
            # chromosome, extent of the chromosome
            yield seqid, ranks[0][1], ranks[-1][1]

class BedpeLine(object):
    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid1 = args[0]
        self.start1 = int(args[1]) + 1
        self.end1 = int(args[2])
        self.seqid2 = args[3]
        self.start2 = int(args[4]) + 1
        self.end2 = int(args[5])
        self.accn = args[6]
        self.score = args[7]
        self.strand1 = args[8]
        self.strand2 = args[9]
        self.isdup = False

    @property
    def innerdist(self):
        if self.seqid1 != self.seqid2:
            return -1
        return abs(self.start2 - self.end1)

    @property
    def outerdist(self):
        if self.seqid1 != self.seqid2:
            return -1
        return abs(self.end2 - self.start1)

    @property
    def is_innie(self):
        return (self.strand1, self.strand2) == ('+', '-')

    def rc(self):
        self.strand1 = '+' if self.strand1 == '-' else '-'
        self.strand2 = '+' if self.strand2 == '-' else '-'

    def _extend(self, rlen, size, start, end, strand):
        if strand == '+':
            end = start + rlen - 1
            if end > size:
                end = size
                start = end - rlen + 1
        else:
            start = end - rlen + 1
            if start < 1:
                start = 1
                end = start + rlen - 1
        return start, end, strand

    def extend(self, rlen, size):
        self.start1, self.end1, self.strand1 = self._extend(\
                rlen, size, self.start1, self.end1, self.strand1)
        self.start2, self.end2, self.strand2 = self._extend(\
                rlen, size, self.start2, self.end2, self.strand2)

    def __str__(self):
        args = (self.seqid1, self.start1 - 1, self.end1,
                self.seqid2, self.start2 - 1, self.end2,
                self.accn, self.score, self.strand1, self.strand2)
        return "\t".join(str(x) for x in args)

    @property
    def bedline(self):
        assert self.seqid1 == self.seqid2
        assert self.start1 <= self.end2
        args = (self.seqid1, self.start1 - 1, self.end2, self.accn)
        return "\t".join(str(x) for x in args)

class BedEvaluate (object):
    def __init__(self, TPbed, FPbed, FNbed, TNbed):

        self.TP = Bed(TPbed).sum(unique=True)
        self.FP = Bed(FPbed).sum(unique=True)
        self.FN = Bed(FNbed).sum(unique=True)
        self.TN = Bed(TNbed).sum(unique=True)

    def __str__(self):
        from jcvi.utils.table import tabulate

        table = {}
        table[("Prediction-True", "Reality-True")] = self.TP
        table[("Prediction-True", "Reality-False")] = self.FP
        table[("Prediction-False", "Reality-True")] = self.FN
        table[("Prediction-False", "Reality-False")] = self.TN
        msg = str(tabulate(table))

        msg += "\nSensitivity [TP / (TP + FN)]: {0:.1f} %\n".\
                format(self.sensitivity * 100)
        msg += "Specificity [TP / (TP + FP)]: {0:.1f} %\n".\
                format(self.specificity * 100)
        msg += "Accuracy [(TP + TN) / (TP + FP + FN + TN)]: {0:.1f} %".\
                format(self.accuracy * 100)
        return msg

    @property
    def sensitivity(self):
        if self.TP + self.FN == 0:
            return 0
        return self.TP * 1. / (self.TP + self.FN)

    @property
    def specificity(self):
        if self.TP + self.FP == 0:
            return 0
        return self.TP * 1. / (self.TP + self.FP)

    @property
    def accuracy(self):
        if self.TP + self.FP + self.FN + self.TN == 0:
            return 0
        return (self.TP + self.TN) * 1. / \
               (self.TP + self.FP + self.FN + self.TN)

    @property
    def score(self):
        return "|".join(("{0:.3f}".format(x) for x in \
                    (self.sensitivity, self.specificity, self.accuracy)))

class BedSummary(object):
    def __init__(self, bed):
        mspans = [(x.span, x.accn) for x in bed]
        spans, accns = zip(*mspans)
        self.mspans = mspans
        self.stats = SummaryStats(spans)
        self.nseqids = len(set(x.seqid for x in bed))
        self.nfeats = len(bed)
        self.total_bases = bed_sum(bed, unique=False)
        self.unique_bases = bed_sum(bed)
        self.coverage = self.total_bases * 1. / self.unique_bases

    def report(self):
        print >> sys.stderr, "Total seqids: {0}".format(self.nseqids)
        print >> sys.stderr, "Total ranges: {0}".format(self.nfeats)
        print >> sys.stderr, "Total unique bases: {0} bp".format(thousands(self.unique_bases))
        print >> sys.stderr, "Total bases: {0} bp".format(thousands(self.total_bases))
        print >> sys.stderr, "Estimated coverage: {0:.1f}x".format(self.coverage)
        print >> sys.stderr, self.stats
        maxspan, maxaccn = max(self.mspans)
        minspan, minaccn = min(self.mspans)
        print >> sys.stderr, "Longest: {0} ({1})".format(maxaccn, maxspan)
        print >> sys.stderr, "Shortest: {0} ({1})".format(minaccn, minspan)

    def __str__(self):
        return "\t".join(str(x) for x in (self.nfeats, self.unique_bases))

def size(args):
    size = 0
    for b in Bed(args.fi):
        size += b.end - b.start + 1
    print(size)

def bed_sum(beds, seqid=None, unique=True):
    if seqid:
        ranges = [(x.seqid, x.start, x.end) for x in beds \
                    if x.seqid == seqid]
    else:
        ranges = [(x.seqid, x.start, x.end) for x in beds]

    unique_sum = range_union(ranges)
    raw_sum = sum(x.span for x in beds)
    return unique_sum if unique else raw_sum

def flanking(args):
    """
    %prog flanking bedfile [options]

    Get up to n features (upstream or downstream or both) flanking a given position.
    """
    from numpy import array, argsort

    p = OptionParser(flanking.__doc__)
    p.add_option("--chrom", default=None, type="string",
            help="chrom name of the position in query. Make sure it matches bedfile.")
    p.add_option("--coord", default=None, type="int",
            help="coordinate of the position in query.")
    p.add_option("-n", default=10, type="int",
            help="number of flanking features to get [default: %default]")
    p.add_option("--side", default="both", choices=("upstream", "downstream", "both"),
            help="which side to get flanking features [default: %default]")
    p.add_option("--max_d", default=None, type="int",
            help="features <= max_d away from position [default: %default]")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if any([len(args) != 1, opts.chrom is None, opts.coord is None]):
        sys.exit(not p.print_help())

    bedfile, = args
    position = (opts.chrom, opts.coord)
    n, side, maxd = opts.n, opts.side, opts.max_d

    chrombed = Bed(bedfile).sub_bed(position[0])

    if side == "upstream":
        data = [(abs(f.start-position[1]), f) for f in chrombed \
            if f.start <= position[1]]
    elif side == "downstream":
        data = [(abs(f.start-position[1]), f) for f in chrombed \
            if f.start >= position[1]]
    else:
        data = [(abs(f.start-position[1]), f) for f in chrombed]

    if maxd:
        data = [f for f in data if f[0]<=maxd]

    n += 1 # not counting self
    n = min(n, len(data))
    distances, subbed = zip(*data)
    distances = array(distances)
    idx = argsort(distances)[:n]
    flankingbed = [f for (i, f) in enumerate(subbed) if i in idx]

    fw = must_open(opts.outfile, "w")
    for atom in flankingbed:
        print >>fw, str(atom)

    return (position, flankingbed)

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
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


    
    
"""('depth', 'calculate average depth per feature using coverageBed'),
        ('mergebydepth', 'returns union of features beyond certain depth'),
        ('sort', 'sort bed file'),
        ('merge', 'merge bed files'),
        ('index', 'index bed file using tabix'),
        ('bins', 'bin bed lengths into each window'),
        ('summary', 'summarize the lengths of the intervals'),
        ('evaluate', 'make truth table and calculate sensitivity and specificity'),
        ('pile', 'find the ids that intersect'),
        ('pairs', 'estimate insert size between paired reads from bedfile'),
        ('mates', 'print paired reads from bedfile'),
        ('sizes', 'infer the sizes for each seqid'),
        ('uniq', 'remove overlapping features with higher scores'),
        ('longest', 'select longest feature within overlapping piles'),
        ('bedpe', 'convert to bedpe format'),
        ('distance', 'calculate distance between bed features'),
        ('sample', 'sample bed file and remove high-coverage regions'),
        ('refine', 'refine bed file using a second bed file'),
        ('flanking', 'get n flanking features for a given position'),
        ('some', 'get a subset of bed features given a list'),
        ('fix', 'fix non-standard bed files'),
        ('filter', 'filter bedfile to retain records between size range'),
        ('filterbedgraph', 'filter bedgraph to extract unique regions'),
        ('random', 'extract a random subset of features'),
        ('juncs', 'trim junctions.bed overhang to get intron, merge multiple beds'),
        ('seqids', 'print out all seqids on one line'),
        ('alignextend', 'alignextend based on BEDPE and FASTA ref'),
        ('clr', 'extract clear range based on BEDPE'),
        ('chain', 'chain bed segments together'),
        ('density', 'calculates density of features per seqid'),
        ('tiling', 'compute the minimum tiling path')
"""

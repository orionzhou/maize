#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column
import pysam
import HTSeq
import collections

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob, create_job_chain

class BamStat(object):
    def __init__(self):
        stats = '''
            pair unpair
            pair_bad pair_dup
            unpair_bad unpair_dup
            pair_map pair_orphan pair_unmap
            unpair_map unpair_unmap 
            pair_map_hq pair_orphan_hq unpair_map_hq
            pair_map0 pair_orphan0 unpair_map0
            pair_map_hq0 pair_orphan_hq0 unpair_map_hq0'''.split()
        for stat in stats:
            setattr(self, stat, 0)
        self.stats = stats
        self.idic = dict()
        self.rdic = dict()

    def __str__(self):
        return "\n".join(["%s\t%d" % (s, getattr(self, s)) for s in self.stats])+"\n"

    __repr__ = __str__

    def write(self, outfile):
        fh = must_open(outfile, "w")
        fh.write(str(self))
        fh.close()

def htseq_read_gtf(fg):
    gtf = HTSeq.GFF_Reader(fg)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded = True)
    for feat in gtf:
        if feat.type == 'exon':
            exons[feat.iv] += feat.attr['gene_id']

def bam_count(args):
    bam = HTSeq.SAM_Reader(args.fi)
    #exons = htseq_read_gtf(args.fg)
    cnts = collections.Counter()
    for bundle in HTSeq.pair_SAM_alignments_with_buffer(bam):
        if len(bundle) != 1:
            continue
        aln1, aln2 = bundle[0]
        if not aln1.aligned and aln2.aligned:
            cnts["_unmapped"] += 1
            continue
        gids = set()
        for iv, val in exons[aln1.iv].steps():
            gids |= val
        for iv, val in exons[aln2.iv].steps():
            gids |= val
        if len(gids) == 1:
            gid = list(gids)[0]
            cnts[gid] += 1
        elif len(gids) == 0:
            cnts["_no_feature"] += 1
        else:
            cnts["_ambiguous"] += 1
    for gid in cnts:
        print("%s\t%d" % (gid, cnts[gid]))

if __name__ == "__main__":
 
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'count feature-overlapping reads in BAM'
    )
    parser.add_argument('fi', help = 'input SAM/BAM file')
    parser.add_argument('fg', help = 'gene GTF file')
    args = parser.parse_args()
    bam_count(args)

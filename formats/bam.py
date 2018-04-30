#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column
import pysam

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

def check_bam(fbam):
    exist_bam = 1
    try:
        bam = pysam.AlignmentFile(fbam, "rb")
    except:
        exist_bam = 0
    return exist_bam

def is_perfect_match(aln):
    nts, bks = aln.get_cigar_stats()
    dst = nts[1]+nts[2]+nts[8]+nts[10] # IDX NM
    if dst == 0:
        return True
    else:
        return False
    
def bam_stat(args):
    bam = pysam.AlignmentFile(args.fi)
    s = BamStat()
    for aln in bam:
        if aln.is_secondary or aln.is_supplementary:
            continue
        if aln.is_paired:
            if aln.is_read1:
                s.pair += 1
            if aln.is_qcfail:
                if aln.is_read1:
                    s.pair_bad += 1
                continue
            if aln.is_duplicate:
                if aln.is_read1:
                    s.pair_dup += 1
                continue
            if aln.is_unmapped and aln.mate_is_unmapped:
                if aln.is_read1:
                    s.pair_unmap += 1
            elif not aln.is_unmapped and not aln.mate_is_unmapped:
                if aln.is_read1:
                    s.pair_map += 1
                aid = aln.query_name
                mq = aln.mapping_quality
                pm = is_perfect_match(aln)
                if aid not in s.rdic:
                    s.rdic[aid] = [mq, pm]
                    continue
                else:
                    m_mq, m_pm = s.rdic[aid]
                    if pm and m_pm:
                        s.pair_map0 += 1
                    if mq >= 20 or m_mq >= 20:
                        s.pair_map_hq += 1
                        if pm and m_pm:
                            s.pair_map_hq0 += 1
                        if aln.is_proper_pair:
                            isize = abs(aln.template_length)
                            if isize not in s.idic:
                                s.idic[isize] = 0
                            s.idic[isize] += 1
                    del s.rdic[aid]
            elif not aln.is_unmapped:
                s.pair_orphan += 1
                mq = aln.mapping_quality
                pm = is_perfect_match(aln)
                if pm:
                    s.pair_orphan0 += 1
                if mq >= 20:
                    s.pair_orphan_hq += 1
                    if pm:
                        s.pair_orphan_hq0 += 1
            else:
                assert aln.is_unmapped and not aln.mate_is_unmapped, "error99"
        else:
            s.unpair += 1
            if aln.is_qcfail:
                s.unpair_bad += 1
                continue
            if aln.is_duplicate:
                s.unpair_dup += 1
                continue
            if aln.is_unmapped:
                s.unpair_unmap += 1
            else:
                s.unpair_map += 1
                mq = aln.mapping_quality
                pm = is_perfect_match(aln)
                if pm:
                    s.unpair_map0 += 1
                if mq >= 20:
                    s.unpair_map_hq += 1
                    if pm:
                        s.unpair_map_hq0 += 1
    if len(s.rdic) > 0:
        logging.debug("%d 'paired' reads don't have a mate" % len(s.rdic))
    if args.isize:
        fho = must_open(args.isize, "w")
        fho.write("insert_size\tcount\n")
        for ins, cnt in s.idic.items():
            fho.write("%d\t%d\n" % (ins, cnt))
        fho.close()
    for stat in s.stats:
        print("%s\t%d" % (stat, getattr(s, stat)))
 
def bam_binstat(args):
    bam = pysam.AlignmentFile(args.fi)
    cdic = dict()
    for aln in bam:
        if aln.reference_id == -1:
            seqid = 'unmapped'
            beg = 1
        else:
            seqid = aln.reference_name
            beg = aln.reference_start
        mq = aln.mapping_quality
        if seqid not in cdic:
            cdic[seqid] = dict()
        pbin = int(beg / 1000000) + 1
        if pbin not in cdic[seqid]:
            cdic[seqid][pbin] = [0, 0]
        if mq >= 20:
            cdic[seqid][pbin][0] += 1
        else:
            cdic[seqid][pbin][1] += 1
    bam.close()

    for seqid, idic in sorted(cdic.items()):
        for pbin, cnts in sorted(idic.items()):
            cntstr = "\t".join([str(x) for x in cnts])
            print("%s\t%d\t%s" % (seqid, pbin, cntstr))

if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'BAM utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("stat",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'BAM statistics'
    )
    sp1.add_argument('fi', help = 'input SAM/BAM file')
    sp1.add_argument('--isize', help = 'output insert size distribution to file')
    sp1.set_defaults(func = bam_stat)

    sp1 = sp.add_parser("binstat",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'count reads in bins'
    )
    sp1.add_argument('fi', help = 'input SAM/BAM file')
    sp1.add_argument('--ws', default = 1, help = 'window size - 1: 1Mb, 2: 100Kb')
    sp1.add_argument('--mq', default = 1, help = 'mapping quality threshold- 1: 20, 2: 10')
    sp1.set_defaults(func = bam_binstat)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


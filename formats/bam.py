#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import numpy as np
from string import Template
import pysam

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

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
        return "\n".join(["\t".join((s, str(getattr(self, s)))) for s in self.stats])

    __repr__ = __str__

    def write(self, outfile):
        fh = must_open(outfile, "w")
        print(self, file=fh)
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
    dst = nts[1]+nts[2]+nts[8] # IDX
    mm = 0
    if aln.has_tag("NM"):
        mm = aln.get_tag("NM")
    elif aln.has_tag("nM"): # STAR
        mm = aln.get_tag("nM")
    else:
        logging.log("no NM/nM tag")
        sys.exit(1)
    dst += mm
    if dst == 0:
        return True
    else:
        return False

def bam_stat(args):
    bam = pysam.AlignmentFile(args.fi, 'r')

    if not args.bychr:
        s = BamStat()
        for aln in bam:
            count_read(aln, s)

        if len(s.rdic) > 0:
            logging.debug("%d 'paired' reads don't have a mate" % len(s.rdic))

        if args.isize:
            fho = must_open(args.isize, "w")
            print("\t".join(('insert_size','count')), file=fho)
            for ins, cnt in s.idic.items():
                print("%d\t%d\n" % (ins, cnt), file=fho)
            fho.close()

        print(s)
    else:
        ss = dict()
        for ist in bam.get_index_statistics():
            ss[ist.contig] = BamStat()
        for aln in bam:
            if aln.is_unmapped: continue
            chrom = aln.reference_name
            count_read(aln, ss[chrom])
        for chrom, s in ss.items():
            for k in s.stats:
                print("\t".join((chrom, k, str(getattr(s, k)))))


def count_read(aln, s):
    if aln.is_secondary or aln.is_supplementary:
        return
    if aln.is_paired:
        if aln.is_read1:
            s.pair += 1
        if aln.is_qcfail:
            if aln.is_read1:
                s.pair_bad += 1
            return
        if aln.is_duplicate:
            if aln.is_read1:
                s.pair_dup += 1
            return
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
                return
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
            return
        if aln.is_duplicate:
            s.unpair_dup += 1
            return
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
    return

def bam_filter(args):
    ibam = pysam.AlignmentFile(args.fi, "rb")
    obam = pysam.AlignmentFile(args.fo, "wb", template=ibam)
    cu = 0
    for aln in ibam:
        if aln.is_secondary or aln.is_supplementary:
            continue
        if aln.is_paired:
            if aln.is_qcfail:
                continue
            if aln.is_duplicate:
                continue
            if aln.is_unmapped:
                continue
        else:
            if aln.is_qcfail:
                continue
            if aln.is_duplicate:
                continue
            if aln.is_unmapped:
                continue
        aid = aln.query_name
        mq = aln.mapping_quality
        pm = is_perfect_match(aln)
        if not pm and mq >= 20:
            obam.write(aln)
            cu += 1
    print("%d hq reads with mismatches" % cu) 

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
    sp1.add_argument('--bychr', action='store_true', help='report stats for each chrom')
    sp1.set_defaults(func = bam_stat)

    sp1 = sp.add_parser("filter",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'filter SAM/BAM file'
    )
    sp1.add_argument('fi', help='input SAM/BAM file')
    sp1.add_argument('fo', help='output BAM file')
    sp1.set_defaults(func = bam_filter)

    sp1 = sp.add_parser("binstat",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'count reads in bins'
    )
    sp1.add_argument('fi', help='input SAM/BAM file')
    sp1.add_argument('--ws', default=1, help='window size - 1: 1Mb, 2: 100Kb')
    sp1.add_argument('--mq', default=1, help='mapping quality threshold- 1: 20, 2: 10')
    sp1.set_defaults(func = bam_binstat)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import pysam

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

def sam2tsv(args):
    sMatch, sMisMatch, sGapOpen, sGapExtend = 2, -3, -5, -2
    sam = pysam.AlignmentFile(args.fi, "r")
   
    print("qId\tqBeg\tqEnd\tqSrd\tqSize\ttId\ttBeg\ttEnd\ttSrd\ttSize\t" +
            "alnLen\tmatch\tmisMatch\tbaseN\tqNumIns\ttNumIns\tqBaseIns\ttBaseIns\tident\tscore\t" +
            "qLoc\ttLoc")
    for x in sam.fetch():
        if x.is_unmapped:
            continue
        tId, tBeg, tEnd, tSrd, tSize = x.reference_name, x.reference_start, x.reference_end, "+", x.reference_length
        qId, qBeg, qEnd, qSrd, qSize = x.query_name, x.query_alignment_start, x.query_alignment_end, "+", x.query_length
        tBeg += 1
        qBeg += 1
        if x.is_reverse: qSrd = "-"
        if args.paired:
            if x.is_read2:
                qId += ".2"
            else:
                qId += ".1"
        alnLen, match, misMatch, baseN, qLen = 0,0,0,0,0
        qNumIns, tNumIns, qBaseIns, tBaseIns = 0,0,0,0
        for op, nt in x.cigartuples:
            if op == 0 or op == 7 or op == 8: # M=X
                alnLen += nt
                qLen += nt
            elif op == 1: # I
                qNumIns += 1
                qBaseIns += nt
                qLen += nt
            elif op == 4: # S
                qLen += nt
            elif op == 2 or op == 3: # DN
                tNumIns += 1
                tBaseIns += nt
            elif op == 5:
                logging.error("hard clipping: %s -> %s:%d" % (qId, tId, tBeg))
                sys.exit(1)
        if x.has_tag("NM"):
            misMatch = x.get_tag("NM")
        match = alnLen - misMatch
        if qSize == 0:
            qSize = qLen
        #assert qSize == qEnd, "error qSize: %d > %d" % (qSize, qEnd)
        
        score_match = match * sMatch
        score_misMatch = misMatch * sMisMatch
        numIns = qNumIns + tNumIns
        score_indel = 0
        if numIns >= 1:
            score_indel = sGapOpen + (numIns - 1) * sGapExtend
        score = score_match + score_misMatch + score_indel
        ident = match / (match + misMatch)

        print("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.03f\t%d\t%s\t%s" % 
                (qId, qBeg, qEnd, qSrd, qSize,
                 tId, tBeg, tEnd, tSrd, tSize,
                 alnLen, match, misMatch, 0,
                 qNumIns, tNumIns, qBaseIns, tBaseIns, ident, score, '', ''))
        if 1 < 2:
            sys.exit(1)

def count(args):
    from math import ceil
    winsize, min_qual = args.winsize, args.min_qual
    cdic = dict()
    bam = pysam.AlignmentFile(args.fi, "rb")
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
        pbin = ceil(beg / winsize)
        if pbin not in cdic[seqid]:
            cdic[seqid][pbin] = [0, 0]
        if mq >= min_qual:
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
            description = 'sam utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv", help = "sam -> tsv")
    sp1.add_argument('fi', help = 'input *.sam or *.bam file')
    sp1.add_argument('--paired', action = "store_true", help = 'paired end input ')
    sp1.set_defaults(func = sam2tsv)
 
    sp1 = sp.add_parser("count",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'count reads in bins'
    )
    sp1.add_argument('fi', help='input *.sam or *.bam file')
    sp1.add_argument('--winsize', default=100000, choices=[100000, 1000000], help='window size')
    sp1.add_argument('--min_qual', default=20, choices=[0, 10, 20], help='min read mapping quality')
    sp1.set_defaults(func = count)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



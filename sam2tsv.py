#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import argparse
import pysam

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'Convert SAM file to TSV'
    )
    parser.add_argument(
            'fi', help = 'SAM/BAM file'
    )
    parser.add_argument(
            'fo', help = 'output file (tsv)'
    )
    args = parser.parse_args()
    fi, fo = args.fi, args.fo
    sMatch, sMisMatch, sGapOpen, sGapExtend = 2, -3, -5, -2
    sam = pysam.AlignmentFile(fi, "r")
   
    fho = open(fo, "w")
    fho.write("qId\tqBeg\tqEnd\tqSrd\tqSize\ttId\ttBeg\ttEnd\ttSrd\ttSize\t" +
            "alnLen\tmatch\tmisMatch\tbaseN\tqNumIns\ttNumIns\tqBaseIns\ttBaseIns\tident\tscore\t" +
            "qLoc\ttLoc\n")
    for x in sam.fetch():
        if x.is_unmapped:
            continue
        tId, tBeg, tEnd, tSrd, tSize = x.reference_name, x.reference_start, x.reference_end, "+", x.reference_length
        qId, qBeg, qEnd, qSrd, qSize = x.query_name, x.query_alignment_start, x.query_alignment_end, "+", x.query_length
        tBeg += 1
        qBeg += 1
        if x.is_reverse: qSrd = "-"
        if x.is_read2:
            qId += ".1"
        else:
            qId += ".2"
        alnLen, match, misMatch, baseN = 0,0,0,0
        qNumIns, tNumIns, qBaseIns, tBaseIns = 0,0,0,0
        for op, nt in x.cigartuples:
            if op == 0 or op == 7 or op == 8:
                alnLen += nt
            elif op == 1:
                qNumIns += 1
                qBaseIns += nt
            elif op == 2 or op == 3:
                tNumIns += 1
                tBaseIns += nt
        if x.has_tag("NM"):
            misMatch = x.get_tag("NM")
        match = alnLen - misMatch
        
        score_match = match * sMatch
        score_misMatch = misMatch * sMisMatch
        numIns = qNumIns + tNumIns
        score_indel = 0
        if numIns >= 1:
            score_indel = sGapOpen + (numIns - 1) * sGapExtend
        score = score_match + score_misMatch + score_indel
        ident = match / (match + misMatch)

        fho.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.03f\t%d\t%s\t%s\n" % 
                (qId, qBeg, qEnd, qSrd, qSize,
                 tId, tBeg, tEnd, tSrd, tSize,
                 alnLen, match, misMatch, 0,
                 qNumIns, tNumIns, qBaseIns, tBaseIns, ident, score, '', ''))
        #exit(1)
    fho.close()
    

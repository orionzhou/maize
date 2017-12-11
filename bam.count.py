#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import argparse
import configparser
from string import Template
import pysam

def check_bam(fbam):
    exist_bam = 1
    try:
        bam = pysam.AlignmentFile(fbam, "rb")
    except:
        exist_bam = 0
    return exist_bam

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'Count reads in bins'
    )
    parser.add_argument(
            'fb', help = 'BAM file'
    )
    parser.add_argument(
            'fo', help = 'output file (tsv)'
    )
    parser.add_argument(
            '--ws', default = 1, 
            help = 'window size - 1: 1Mb, 2: 100Kb (default: 1)'
    )
    parser.add_argument(
            '--mq', default = 1, 
            help = 'mapping quality threshold- 1: 20, 2: 10 (default: 1)'
    )
    args = parser.parse_args()

    cdic = dict()
    bam = pysam.AlignmentFile(args.fb, "rb")
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

    fho = open(args.fo, "w")
    for seqid, idic in sorted(cdic.items()):
        for pbin, cnts in sorted(idic.items()):
            cntstr = "\t".join([str(x) for x in cnts])
            fho.write("%s\t%d\t%s\n" % (seqid, pbin, cntstr))
    fho.close()

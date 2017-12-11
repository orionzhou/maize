#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import argparse
import configparser
from string import Template
from colorama import init, Fore, Back, Style
import pysam

def get_read_base(aln, qseq, qual, pos):
    for qpos, rpos, rbase in aln:
        if rpos != pos:
            continue
        return [rbase, qseq[qpos], qual[qpos]]
    return ['', '', 0]
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'Allele-specific Expression calculation'
    )
    parser.add_argument(
            'fb', help = 'bam file'
    )
    parser.add_argument(
            'fo', help = 'output file (tsv)'
    )
    parser.add_argument(
            'fv', nargs = '?', default = '/home/springer/zhoux379/data/misc2/mo17vnt/53.vnt.final/61.rna.bed', help = 'variant file'
    )
    args = parser.parse_args()
    fb, fo, fv = args.fb, args.fo, args.fv
    bam = pysam.AlignmentFile(fb, "rb")

    fhv = open(fv, "r")
    fho = open(fo, "w")
    vdic = dict()
    rset = set()
    for line in fhv:
        line = line.strip("\n")
        seqid, beg, end, gts = line.split("\t")
        beg, end = int(beg), int(end)
        ref, alt = gts.split(",")
        for x in bam.fetch(seqid, beg, end):
            rid = x.query_name
            #if rid in rset: continue
            #if rid != 'D00635:197:CAC47ANXX:2:2209:12819:101106':
            #    continue
            aln = x.get_aligned_pairs(matches_only = True, with_seq = True)
            rbase, qbase, qual = get_read_base(aln, x.query_sequence, x.query_qualities, beg)
            rset.add(rid)
            pair = 'pe' if x.is_paired else 'se'
            first = 'r1' if x.is_read1 else 'r2'
            dup = 'dup' if x.is_secondary else 'uniq'
            fho.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\n" % 
                    (seqid, beg, end, gts, rid, rbase, qbase, qual))
            #exit(1)
        pos = "%s_%d" % (seqid, int(beg) + 1)
        vdic[pos] = [ref, alt]
    fhv.close()
    fho.close()
    

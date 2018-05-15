#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
from string import Template
import pysam

def get_read_base(aln, qseq, qual, pos):
    for qpos, rpos, rbase in aln:
        if rpos != pos:
            continue
        return [rbase, qseq[qpos], qual[qpos]]
    return ['', '', 0]
def infer_variant(x, minBaseQual):
    aln = x.get_aligned_pairs(matches_only = True, with_seq = True)
    qseq = x.query_sequence
    qual = x.query_qualities
    mms = dict()
    for qpos, rpos, rbase in aln:
        if rbase.islower() and qual[qpos] >= minBaseQual:
            mms[qpos] = [rpos+1, 'M', qseq[qpos]]
    res = []
    qpos, rpos = 0, x.reference_start
    rbeg = rpos
    vnts = []
    for opt, nbase in x.cigartuples:
        if opt == 0: #M
            for i in range(0, nbase):
                if qpos+i in mms:
                    vnts.append(mms[qpos+i])
            qpos += nbase
            rpos += nbase
        elif opt == 1: #I
            vnts.append([rpos+1, "I", qseq[qpos+1:qpos+nbase+1]])
            qpos += nbase
        elif opt == 2: #D
            vnts.append([rpos+1, "D", nbase])
            rpos += nbase
        elif opt == 3: #N
            res.append([rbeg, rpos, vnts])
            rpos += nbase
            rbeg = rpos
            vnts = []
        elif opt == 4: #S
            qpos += nbase
        elif opt == 7 or opt == 8: #=X
            qpos += 1
            rpos += 1
    if rpos > rbeg:
        res.append([rbeg, rpos, vnts])
    if qpos != x.query_length or rpos != x.reference_end:
        print(qpos, x.query_alignment_end, rpos, x.reference_end, x.cigartuples, aln)
        exit(1)
    return res
def read_variants(fv):
    fhv = open(fv, "r")
    vdic = dict()
    for line in fhv:
        line = line.strip("\n")
        seqid, beg, end, gts = line.split("\t")
        beg, end = int(beg), int(end)
        ref, alt = gts.split(",")
        pos = "%s_%d" % (seqid, int(beg) + 1)
        vdic[pos] = [ref, alt]
    fhv.close()
    return vdic
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            description = 'Allele-specific Expression calculation'
    )
    parser.add_argument(
            'fi', help = 'bam file'
    )
    parser.add_argument(
            'fo', help = 'output file (bed)'
    )
    parser.add_argument(
            '--min_mapq', default = 20, help = 'min mapping quality [20]'
    )
    parser.add_argument(
            '--min_baseq', default = 20, help = 'min base quality [20]'
    )
    parser.add_argument(
            '--vcf', nargs = '?', default = '/home/springer/zhoux379/data/misc2/mo17vnt/53.vnt.final/61.rna.bed', help = 'variant file'
    )
    args = parser.parse_args()
    fi, fo = args.fi, args.fo
    min_mapq, min_baseq = args.min_mapq, args.min_baseq
    bam = pysam.AlignmentFile(fi, "rb")
    #vnts = pybedtools.BedTool(args.vcf)
    
    fho = open(fo, "w")
    for x in bam.fetch():
        if x.is_duplicate or x.is_unmapped or x.is_secondary or x.mapping_quality < min_mapq:
            continue
        sid = x.query_name
        rid, rbeg, rend = x.reference_name, x.reference_start, x.reference_end
        #pair = 'pe' if x.is_paired else 'se'
        #first = 'r1' if x.is_read1 else 'r2'
        res = infer_variant(x, min_baseq)
        #if sid == 'D00635:197:CAC47ANXX:2:2311:18289:97712':
        #    print(res)
        #    exit(1)
        for rbeg, rend, vnts in res: 
            vntstr = '.'
            if len(vnts) > 0:
                vntstr = " ".join([":".join(map(str,vnt)) for vnt in vnts])
            fho.write("%s\t%d\t%d\t%s\t%s\n" % 
                    (rid, rbeg, rend, sid, vntstr))
        #exit(1)
    fho.close()
    

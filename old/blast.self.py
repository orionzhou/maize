#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math
import os.path as op
import numpy as np
import argparse

def blast2tbl(fi, fo):
    e_min = 1
    fhi = open(fi, "r")
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, qlen, tid, tlen, alen, bit, e) = line.split("\t")
        e = float(e)
        if e != 0 and e < e_min: e_min = e
    fhi.close()

    fhi = open(fi, "r")
    fho = open(fo, "w")
    (cqid, ctid, cqlen, ctlen, calen, cscore) = ('', '', 0, 0, 0, 0)
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, qlen, tid, tlen, alen, bit, e) = line.split("\t")
        if qid == tid: continue
        (qlen, tlen, alen) = (float(qlen), float(tlen), float(alen))
        (e, bit) = (float(e), float(bit))
        if e == 0: e = e_min
        score = - math.log10(e)
        if qid == cqid and tid == ctid:
            calen += alen
            cscore += score / 2
        else:
            if cqid != '':
                if calen/min(cqlen, ctlen) >= 0.5 and cscore > 10:
                    print >>fho, "%s\t%s\t%g" % (cqid, ctid, cscore)
            (cqid, ctid, cqlen, ctlen, calen, cscore) = (qid, tid, qlen, tlen, alen, score)
    if calen/min(cqlen, ctlen) >= 0.5 and cscore > 10:
        print >>fho, "%s\t%s\t%g" % (cqid, ctid, cscore)
    fhi.close()
    fho.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "blast against itself"
    )
    parser.add_argument(
        'fi', help = "input file (fasta)"
    )
    parser.add_argument(
        'fo', help = "output file (prefix)"
    )
    nproc = int(os.environ['nproc'])
    parser.add_argument(
        '--cpu', dest="cpu", default=nproc, type=int,
        help='number processors to use (default: all/%d)' % nproc)
    args = parser.parse_args()
    (fi, fo) = (args.fi, args.fo)
    ncpu = args.cpu
    
    cmd = "makeblastdb -dbtype prot -in %s -out %s.1" % (fi, fo)
    os.system(cmd)
    f01 = fo + ".1.blast.tbl"
    cmd = "blastp -db %s.1 -query %s -out %s.out -outfmt 6 -max_target_seqs 1000 -evalue 1e-5 -num_threads %d" % (fo, fi, fo, ncpu)
    print cmd
    os.system(cmd)


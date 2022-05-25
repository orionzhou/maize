#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MCL clustering
"""
import os
import os.path as op
import sys
import math
import logging

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def blast2tsv(args, fi, fo):
    e_min = 1
    fhi = open(fi, "r")
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, tid, pident, alen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore, qlen, tlen) = line.split("\t")
        e = float(evalue)
        if e != 0 and e < e_min: e_min = e
    fhi.close()

    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, tid, pident, alen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore, qlen, tlen) = line.split("\t")
        (qlen, tlen, alen) = (float(qlen), float(tlen), float(alen))
        (e, bit) = (float(evalue), float(bitscore))
        if e == 0: e = e_min
        rlen = alen / min(qlen, tlen)
        if qid != tid and rlen >= args.min_cov and e < args.max_evalue:
            score = - math.log10(e)
            fho.write("%s\t%s\t%g\n" % (qid, tid, score))
    fhi.close()
    fho.close()

def mcl2tsv(args, fi, fo):
    fhi = open(fi, "r")
    fho = open(fo, "w")
    fho.write("grp\tgid\n")
    grp = 1
    for line in fhi:
        line = line.strip("\n")
        gids = line.split("\t")
        if len(gids) < args.min_cluster_size:
            continue
        for gid in gids:
            fho.write("%d\t%s\n" % (grp, gid))
        grp += 1
    fhi.close()
    fho.close()

def main(args):
    f1 = args.fi
    f2 = f"tmp.1.tsv"
    cmd = f"sort -nk 11 {f1} | sort -uk1,2 > {f2}"
    sh(cmd)

    f3 = "tmp.2.tsv"
    blast2tsv(args, f2, f3)

    f4 = 'tmp.3.txt'
    cmd = f"mcl {f3} -I {args.I} -te {args.thread} --abc -o {f4}"
    sh(cmd)

    f5 = args.fo
    mcl2tsv(args, f4, f5)

    sh(f"rm {f2} {f3} {f4}")

if __name__ == '__main__':
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'sequence clustering using mcl'
    )

    fmtdesc = '''
  1.  qseqid      query or source (e.g., gene) sequence id
  2.  sseqid      subject  or target (e.g., reference genome) sequence id
  3.  pident      percentage of identical matches
  4.  length      alignment length (sequence overlap)
  5.  mismatch    number of mismatches
  6.  gapopen     number of gap openings
  7.  qstart      start of alignment in query
  8.  qend        end of alignment in query
  9.  sstart      start of alignment in subject
 10.  send        end of alignment in subject
 11.  evalue      expect value
 12.  bitscore    bit score
 13.  qlen        query length
 14.  tlen        target length
'''

    ps.add_argument('fi', help=f"blast tabular file (14 columns: \n{fmtdesc})")
    ps.add_argument('fo', help='mcl clustering output')
    ps.add_argument('--type', default="protein", choices=['protein','cds'], help='sequence type')
    ps.add_argument('--min_cov', type=float, default=0.5, help='min hit coverage')
    ps.add_argument('--max_evalue', type=float, default=1e-10, help='max hit evalue')
    ps.add_argument('-I', type=float, default=2.0, help='mcl inflation')
    ps.add_argument('--thread', type=int, default=4, help='mcl threads')
    ps.add_argument('--min_cluster_size', type=int, default=3, help='min cluster size')

    args = ps.parse_args()
    main(args)



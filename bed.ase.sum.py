#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            description = 'Allele-specific Expression calculation'
    )
    parser.add_argument(
            'fr', help = 'input read ASE file (tsv)'
    )
    parser.add_argument(
            'fb', help = 'input gene-read intersection file (bed)'
    )
    parser.add_argument(
            'fo', help = 'output gene ASE file (tsv)'
    )
    args = parser.parse_args()
    fr, fb, fo = args.fr, args.fb, args.fo
    
    fhr = open(fr, "r")
    rdic = dict()
    for line in fhr:
        rid, n0, n1, nunk, nerr = line.strip("\n").split("\t")
        if rid == 'rid':
            continue
        tag = "unk"
        if int(n0) > 0 and int(n1) == 0:
            tag = 'h0'
        elif int(n0) == 0 and int(n1) > 0:
            tag = 'h1'
        elif int(n0) > 0 and int(n1) > 0:
            tag = 'cft'
        rdic[rid] = tag
    fhr.close()

    fhb = open(fb, "r")
    gdic, tdic = dict(), dict()
    for line in fhb:
        row = line.strip("\n").split("\t")
        gid, rid = row[3], row[7]
        if gid not in gdic:
            gdic[gid] = {'n0': 0, 'n1': 0, 'ncft': 0}
        if gid not in tdic:
            tdic[gid] = set()
        if rid not in tdic[gid]:
            tdic[gid].add(rid)
        else:
            continue
        if rid not in rdic:
            print("%s not in read dict" % rid)
        if rdic[rid] == 'h0':
            gdic[gid]['n0'] += 1
        elif rdic[rid] == 'h1':
            gdic[gid]['n1'] += 1
        elif rdic[rid] == 'cft':
            gdic[gid]['ncft'] += 1
    fhb.close()
    
    fho = open(fo, "w")
    fho.write("gid\tn0\tn1\tncft\n")
    for gid in sorted(gdic):
        sdic = gdic[gid]
        n0, n1, ncft = sdic['n0'], sdic['n1'], sdic['ncft']
        if n0 + n1 < 0:
            continue
        fho.write("%s\t%d\t%d\t%d\n" % (gid, n0, n1, ncft))
    fho.close()
    

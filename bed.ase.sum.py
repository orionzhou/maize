#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import argparse
import configparser
from colorama import init, Fore, Back, Style

if __name__ == "__main__":
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
        rid, nref, nalt, nunk, nerr = line.strip("\n").split("\t")
        if rid == 'rid':
            continue
        tag = "unk"
        if int(nref) > 0 and int(nalt) == 0:
            tag = 'ref'
        elif int(nref) == 0 and int(nalt) > 0:
            tag = 'alt'
        elif int(nref) > 0 and int(nalt) > 0:
            tag = 'cft'
        rdic[rid] = tag
    fhr.close()

    fhb = open(fb, "r")
    gdic, tdic = dict(), dict()
    for line in fhb:
        row = line.strip("\n").split("\t")
        gid, rid = row[3], row[7]
        if gid not in gdic:
            gdic[gid] = {'nref': 0, 'nalt': 0, 'ncft': 0}
        if gid not in tdic:
            tdic[gid] = set()
        if rid not in tdic[gid]:
            tdic[gid].add(rid)
        else:
            continue
        if rid not in rdic:
            print("%s not in read dict" % rid)
        if rdic[rid] == 'ref':
            gdic[gid]['nref'] += 1
        elif rdic[rid] == 'alt':
            gdic[gid]['nalt'] += 1
        elif rdic[rid] == 'cft':
            gdic[gid]['ncft'] += 1
    fhb.close()
    
    fho = open(fo, "w")
    fho.write("gid\tnref\tnalt\tncft\n")
    for gid in sorted(gdic):
        sdic = gdic[gid]
        nref, nalt, ncft = sdic['nref'], sdic['nalt'], sdic['ncft']
        if nref + nalt < 0:
            continue
        fho.write("%s\t%d\t%d\t%d\n" % (gid, nref, nalt, ncft))
    fho.close()
    

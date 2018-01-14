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

def infer_read(rows, fho, fhb):
    rdic, vdic = dict(), dict()
    sid, beg, end, rid = rows[0][0:4]
    locs = set()
    for row in rows:
        assert rid == row[3], "error: " + str(rows) 
        if row[4] != '.':
            for rvnt in row[4].split(" "):
                pos, vtype, vnt = rvnt.split(":")
                rdic[pos] = "%s:%s" % (vtype, vnt)
        vpos, vnt = row[7], row[8]
        vdic[vpos] = vnt
        loc = "%s:%s" % tuple(row[0:2])
        if loc not in locs:
            locs.add(loc)
            fhb.write("%s\t%s\t%s\t%s\n" % tuple(row[0:4]))
    nref, nalt, nunk, nerr = 0, 0, 0, 0
    for pos in vdic:
        if pos in rdic:
            if rdic[pos] == vdic[pos]:
                nalt = nalt + 1
            else:
                nunk = nunk + 1
            del rdic[pos]
        else:
            nref = nref + 1
    nerr = len(rdic)
    fho.write("%s\t%d\t%d\t%d\t%d\n" % (rid, nref, nalt, nunk, nerr))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'Allele-specific Expression calculation'
    )
    parser.add_argument(
            'fi', help = 'input bed (ase) file'
    )
    parser.add_argument(
            'fo', help = 'output read ASE file (tsv)'
    )
    parser.add_argument(
            'fb', help = 'output read location file (bed)'
    )
    parser.add_argument(
            '--min_baseq', default = 20, help = 'min base quality [20]'
    )
    args = parser.parse_args()
    fi, fo, fb = args.fi, args.fo, args.fb
    min_baseq = args.min_baseq
    
    fhi = open(fi, "r")
    fho = open(fo, "w")
    fhb = open(fb, "w")
    fho.write("rid\tnref\tnalt\tnunk\tnerr\n")
    prid = ""
    rows = []
    for line in fhi:
        row = line.strip("\n").split("\t")
        if prid == "":
            rows.append(row)
            prid = row[3]
        elif row[3] == prid:
            rows.append(row)
        else:
            infer_read(rows, fho, fhb)
            rows = [row]
            prid = row[3]
    infer_read(rows, fho, fhb)
    fhi.close()
    fho.close()
    fhb.close()
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
from string import Template

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
        vpos, vnt, phase = row[7:10]
        assert phase in ['0|1','1|0'], "Unknown phase: %s" % phase
        vdic[vpos] = [vnt, phase]
        loc = "%s:%s" % tuple(row[0:2])
        if loc not in locs:
            locs.add(loc)
            fhb.write("%s\t%s\t%s\t%s\n" % tuple(row[0:4]))
    n0, n1, nunk, nerr = 0, 0, 0, 0
    for pos in vdic:
        if pos in rdic:
            if rdic[pos] == vdic[pos][0]:
                if vdic[pos][1] == '0|1':
                    n1 = n1 + 1
                else:
                    n0 = n0 + 1
            else:
                nunk = nunk + 1
            del rdic[pos]
        else:
            if vdic[pos][1] == '0|1':
                n0 = n0 + 1
            else:
                n1 = n1 + 1
    nerr = len(rdic)
    fho.write("%s\t%d\t%d\t%d\t%d\n" % (rid, n0, n1, nunk, nerr))
if __name__ == "__main__":
    import argparse
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
    fho.write("rid\tn0\tn1\tnunk\tnerr\n")
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
    

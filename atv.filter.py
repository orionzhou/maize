#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'filter an alignment file'
    )
    parser.add_argument(
        'fi', help = 'input alignment file (tsv)'
    )
    parser.add_argument(
        'fo', help = 'output file (*.tsv)'
    )
    parser.add_argument(
        '--match', type=int, default=1, help='min matches [1]'
    )
    parser.add_argument(
        '--ident', type=float, default=0.8, help='min identity [0.8]'
    )
    parser.add_argument(
        '--cov', type=float, default=0.8, help='min query coverage [0.8]'
    )
    parser.add_argument(
        '--best', action = 'store_true', help='only keep best hit(s) [N]'
    )
    args = parser.parse_args()
    
    if args.fi and args.fi != '-' and args.fi != 'stdin':
        sys.stdin = open(args.fi, "r")
    if args.fo and args.fo != '-' and args.fo != 'stdout':
        sys.stdout = open(args.fo, "w")
   
    line = sys.stdin.readline()
    print(line.strip("\n"))
    pqid = ''
    pscore = 0
    lines = []
    for line in sys.stdin:
        line = line.strip("\n")
        qName, qStart, qEnd, qSize, strand, \
        tName, tStart, tEnd, tSize, \
        alnLen, match, misMatch, baseN, \
        qNumIns, tNumIns, qBaseIns, tBaseIns, \
        ident, score, qLoc, tLoc = line.split("\t")
        #print(line)
        #print(qSize)
        if float(ident) < args.ident: continue
        if int(alnLen)/int(qSize) < args.cov: continue
        if int(match) < args.match: continue
        if pqid == '':
            pqid = qName
            pscore = score
            lines.append(line)
        elif qName != pqid:
            print("\n".join(lines))
            pqid = qName
            pscore = score
            lines = [line]
        else:
            if args.best:
                if score > pscore:
                    lines = [line]
                    pscore = score
                elif score == pscore:
                    lines.append(line)
            else:
                lines.append(line)
    print("\n".join(lines))
 

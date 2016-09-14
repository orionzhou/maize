#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'format conversion'
    )
    parser.add_argument(
        'fi', help = 'input file (*.clstr)'
    )
    parser.add_argument(
        'fo', help = 'output file (*.tbl)'
    )
    args = parser.parse_args()
   
    (fi, fo) = (args.fi, args.fo)
    fhi = open(fi, "r")
    fho = open(fo, "w")
    print >>fho, "grp\tid"

    grp = 0
    p = re.compile('\d+.*>(\S+)\.{3}\s')
    for line in fhi:
        if line[0] == ">":
            grp += 1
            continue
        line = line.strip("\n")
        m = p.match(line)
        print >>fho, "%d\t%s" % (grp, m.group(1))
    fhi.close()
 

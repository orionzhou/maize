#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

from maize.formats.base import BaseFile, read_block, must_open
from maize.apps.base import need_update, which, sh
from maize.utils.location import locAryLen



def rmgap(args):
    firstLine = True
    pid, locs = '', []
    for line in must_open(args.fi):
        line = line.rstrip("\n")
        if not line:
            continue
        ps = line.split()
        assert len(ps) == 12, "not 12 fields: %s" % line
        cid = "\t".join(ps[0:8])
        oStart, oEnd, oSize = int(ps[9]), int(ps[10]), int(ps[11])
        if firstLine:
            pid = cid
            locs.append([oStart, oEnd, oSize])
            firstLine = False
        elif pid == cid:
            locs.append([oStart, oEnd, oSize])
        else:
            rm1gap(pid, locs)
            pid = cid
            locs = [[oStart, oEnd, oSize]]
    rm1gap(pid, locs)

def rm1gap(pid, locs):
    tName, tStart, tEnd, srd, qName, qStart, qEnd, cid = pid.split("\t")
    tStart, tEnd, qStart, qEnd = int(tStart), int(tEnd), int(qStart), int(qEnd)
    if locs[0][2] == 0:
        assert len(locs) == 1 and locs[0][0] == locs[0][1] == -1, "error: %s" % pid
        print(pid)
    else:
        gSize = sum([x[2] for x in locs])
        rlocs = [[max(x[0]-tStart, 0), min(x[1]-tStart, tEnd-tStart)] for x in locs]
        assert locAryLen(rlocs) == gSize, "gap len error: %d %d\n%s" % (gSize, locAryLen(rlocs), pid)
        print("%s\t%d" % (pid, gSize))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'chainBed utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("rmgap", help = "remove genomic gaps from chainBed file")
    sp1.add_argument('fi', help = 'input chainBed file (intersectBed with gapBed)')
    sp1.set_defaults(func = rmgap)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



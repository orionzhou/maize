#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

from jcvi.formats.base import BaseFile, read_block, must_open
from jcvi.formats.chain import ChainLine, Chain
from jcvi.apps.base import need_update, which, sh

def chainstat(args):
    sh("chain.py 2bed %s > tmp.bed" % args.fi)
    logging.debug("total size")
    sh("bed.py size tmp.bed")
    logging.debug("tgt noredundant size")
    sh("cut -f1-3 tmp.bed | sortBed -i stdin | mergeBed -i stdin | bed.py size -")
    logging.debug("qry noredundant size")
    sh("cut -f5-7 tmp.bed | sortBed -i stdin | mergeBed -i stdin | bed.py size -")

def print_chain(cid, tName, qName, qStrand, tSize, qSize, locs, fho):
    chain = "chain"
    score = 1000
    tStrand = "+"
    tStart = min(x[0] for x in locs)
    tEnd = max(x[1] for x in locs)
    qStart = min(x[2] for x in locs)
    qEnd = max(x[3] for x in locs)
    if qStrand == "-":
        qStart, qEnd = qSize - qEnd, qSize - qStart
    headerline = " ".join(str(x) for x in (
         chain, score, tName, tSize, tStrand, tStart,
         tEnd, qName, qSize, qStrand, qStart, qEnd, cid
    ))
    fho.write(headerline + "\n")
    for i in range(len(locs)):
        tb, te, qb, qe = locs[i]
        size, size2 = te - tb, qe - qb
        assert size == size2, "size not equal"
        if i == len(locs) - 1:
            fho.write(f"{size}\n")
        else:
            dt = locs[i+1][0] - te
            if qStrand == "-":
                dq = qb - locs[i+1][3]
            else:
                dq = locs[i+1][2] - qe
            fho.write("%d\t%d\t%d\n" % (size, dt, dq))
    fho.write("\n")

def bed2chain(args):
    from jcvi.formats.sizes import Sizes
    tdic = Sizes(args.tsize)
    qdic = Sizes(args.qsize)

    firstline = True
    cid0, tName0, qName0, srd0, locs = '', '', '', '', []
    fho = must_open(args.fo, 'w')
    for line in must_open(args.fi):
        line = line.rstrip("\n")
        if not line:
            continue
        tName, tStart, tEnd, srd, qName, qStart, qEnd, cid = line.split()[:8]
        tStart, tEnd, qStart, qEnd = int(tStart), int(tEnd), int(qStart), int(qEnd)
        if firstline:
            cid0, tName0, qName0, srd0 = cid, tName, qName, srd
            locs.append([tStart, tEnd, qStart, qEnd])
            firstline = False
        elif cid0 == cid:
            assert tName == tName0 and qName == qName0 and srd == srd0, "inconsistent info in chain"
            locs.append([tStart, tEnd, qStart, qEnd])
        else:
            print_chain(cid0, tName0, qName0, srd0, tdic.get_size(tName0), qdic.get_size(qName0), locs, fho)
            cid0, tName0, qName0, srd0 = cid, tName, qName, srd
            locs = [[tStart, tEnd, qStart, qEnd]]
    print_chain(cid0, tName0, qName0, srd0, tdic.get_size(tName0), qdic.get_size(qName0), locs, fho)

def chain2bed(args):
    chainFile = Chain(args.fi)
    fho = must_open(args.fo, 'w')
    for c in chainFile.chains:
        cid, score, tName, tSize, tSrd, tStart, tEnd, \
                qName, qSize, qSrd, qStart, qEnd, cid = c.chain.split()
        assert tSrd == '+', 'tStrand is not "+"'
        score, tStart, tEnd, tSize, qStart, qEnd, qSize = \
                int(score), int(tStart), int(tEnd), int(tSize), \
                int(qStart), int(qEnd), int(qSize)
        if qSrd == '-':
            qStart, qEnd = qSize - qEnd, qSize - qStart
        offset_t, offset_q = 0, 0
        for ungapped, dt, dq in c.blocks:
            rtb, rte = offset_t, offset_t + ungapped
            rqb, rqe = offset_q, offset_q + ungapped
            tb, te = tStart + rtb, tStart + rte
            if qSrd == '-':
                qb, qe = qEnd - rqe, qEnd - rqb
            else:
                qb, qe = qStart + rqb, qStart + rqe
            tstr = "%s:%d-%d" % (tName, tb, te)
            qstr = "%s:%d-%d" % (qName, qb, qe)
            if args.qry:
                fho.write("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\n" % (qName, qb, qe, qSrd, tName, tb, te, cid))
            else:
                fho.write("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\n" % (tName, tb, te, qSrd, qName, qb, qe, cid))
            offset_t += ungapped + dt
            offset_q += ungapped + dq

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'chain utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2bed", help = "convert to 8-col BED file")
    sp1.add_argument('fi', help = 'input chain file')
    sp1.add_argument('fo', help = 'output bed file')
    sp1.add_argument('--qry', action = 'store_true', help = 'use query coordinate system')
    sp1.set_defaults(func = chain2bed)

    sp1 = sp.add_parser("fromBed", help = "convert from bed to chain file")
    sp1.add_argument('fi', help = 'input bed file')
    sp1.add_argument('tsize', help = 'target size file')
    sp1.add_argument('qsize', help = 'query size file')
    sp1.add_argument('fo', help = 'output chain file')
    sp1.set_defaults(func = bed2chain)

    sp1 = sp.add_parser("stat", help = "get chain stats")
    sp1.add_argument('fi', help = 'input chain file')
    sp1.set_defaults(func = chainstat)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



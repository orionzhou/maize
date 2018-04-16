#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Create the UCSC chain file which is needed to lift over from one coordinate
system to another.

File format:
<http://genome.ucsc.edu/goldenPath/help/chain.html>

chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
  9       1       0
  10      0       5
  48

Header Line:
 chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
Alignment Data Lines
 size dt dq

NOTE: The last line of the alignment section contains only one number: the ungapped
alignment size of the last block.
"""

import os.path as op
import sys
import logging

from maize.formats.base import BaseFile, read_block, must_open
from maize.apps.base import need_update, which

class ChainLine (object):

    def __init__(self, chain, lines):
        self.chain = chain
        self.blocks = []
        for line in lines:
            atoms = line.split()
            if len(atoms) == 1:
                atoms += [0, 0]
            if len(atoms) == 0:
                continue

            self.blocks.append([int(x) for x in atoms])

        self.ungapped, self.dt, self.dq = zip(*self.blocks)
        self.ungapped = sum(self.ungapped)
        self.dt = sum(self.dt)
        self.dq = sum(self.dq)

class Chain (BaseFile):

    def __init__(self, filename):
        super(Chain, self).__init__(filename)
        self.chains = list(self.iter_chain())

        self.ungapped = sum(x.ungapped for x in self.chains)
        self.dt = sum(x.dt for x in self.chains)
        self.dq = sum(x.dq for x in self.chains)

    def __len__(self):
        return len(self.chains)

    def iter_chain(self):
        fp = open(self.filename)
        for chain, lines in read_block(fp, "chain"):
            lines = list(lines)
            yield ChainLine(chain, lines)

def fromagp(args):
    """
    %prog fromagp agpfile componentfasta objectfasta

    Generate chain file from AGP format. The components represent the old
    genome (target) and the objects represent new genome (query).
    """
    from jcvi.formats.agp import AGP
    from jcvi.formats.sizes import Sizes

    p = OptionParser(fromagp.__doc__)
    p.add_option("--novalidate", default=False, action="store_true",
                 help="Do not validate AGP")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    agpfile, componentfasta, objectfasta = args
    chainfile = agpfile.rsplit(".", 1)[0] + ".chain"
    fw = open(chainfile, "w")
    agp = AGP(agpfile, validate=(not opts.novalidate))
    componentsizes = Sizes(componentfasta).mapping
    objectsizes = Sizes(objectfasta).mapping
    chain = "chain"
    score = 1000
    tStrand = "+"
    id = 0
    for a in agp:
        if a.is_gap:
            continue

        tName = a.component_id
        tSize = componentsizes[tName]
        tStart = a.component_beg
        tEnd = a.component_end
        tStart -= 1

        qName = a.object
        qSize = objectsizes[qName]
        qStrand = "-" if a.orientation == "-" else "+"
        qStart = a.object_beg
        qEnd = a.object_end
        if qStrand == '-':
            _qStart = qSize - qEnd + 1
            _qEnd = qSize - qStart + 1
            qStart, qEnd = _qStart, _qEnd
        qStart -= 1

        id += 1
        size = a.object_span
        headerline = "\t".join(str(x) for x in (
             chain, score, tName, tSize, tStrand, tStart,
             tEnd, qName, qSize, qStrand, qStart, qEnd, id
        ))
        alignmentline = size
        print >> fw, headerline
        print >> fw, alignmentline
        print >> fw

    fw.close()
    logging.debug("File written to `{0}`.".format(chainfile))

def faToTwoBit(fastafile):
    twobitfile = fastafile.rsplit(".", 1)[0] + ".2bit"
    cmd = "faToTwoBit {0} {1}".format(fastafile, twobitfile)
    if need_update(fastafile, twobitfile):
        sh(cmd)
    return twobitfile

def blat(args):
    """
    %prog blat old.fasta new.fasta

    Generate psl file using blat.
    """
    p = OptionParser(blat.__doc__)
    p.add_option("--minscore", default=100, type="int",
                 help="Matches minus mismatches gap penalty [default: %default]")
    p.add_option("--minid", default=98, type="int",
                 help="Minimum sequence identity [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    oldfasta, newfasta = args
    twobitfiles = []
    for fastafile in args:
        tbfile = faToTwoBit(fastafile)
        twobitfiles.append(tbfile)

    oldtwobit, newtwobit = twobitfiles
    cmd = "pblat -threads={0}".format(opts.cpus) if which("pblat") else "blat"
    cmd += " {0} {1}".format(oldtwobit, newfasta)
    cmd += " -tileSize=12 -minScore={0} -minIdentity={1} ".\
                format(opts.minscore, opts.minid)
    pslfile = "{0}.{1}.psl".format(*(op.basename(x).split('.')[0] \
                for x in (newfasta, oldfasta)))
    cmd += pslfile
    sh(cmd)

def frompsl(args):
    """
    %prog frompsl old.new.psl old.fasta new.fasta

    Generate chain file from psl file. The pipeline is describe in:
    <http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver>
    """
    from maize.formats.sizes import Sizes

    p = OptionParser(frompsl.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    pslfile, oldfasta, newfasta = args
    pf = oldfasta.split(".")[0]

    # Chain together alignments from using axtChain
    chainfile = pf + ".chain"
    twobitfiles = []
    for fastafile in (oldfasta, newfasta):
        tbfile = faToTwoBit(fastafile)
        twobitfiles.append(tbfile)
    oldtwobit, newtwobit = twobitfiles

    if need_update(pslfile, chainfile):
        cmd = "axtChain -linearGap=medium -psl {0}".format(pslfile)
        cmd += " {0} {1} {2}".format(oldtwobit, newtwobit, chainfile)
        sh(cmd)

    # Sort chain files
    sortedchain = chainfile.rsplit(".", 1)[0] + ".sorted.chain"
    if need_update(chainfile, sortedchain):
        cmd = "chainSort {0} {1}".format(chainfile, sortedchain)
        sh(cmd)

    # Make alignment nets from chains
    netfile = pf + ".net"
    oldsizes = Sizes(oldfasta).filename
    newsizes = Sizes(newfasta).filename
    if need_update((sortedchain, oldsizes, newsizes), netfile):
        cmd = "chainNet {0} {1} {2}".format(sortedchain, oldsizes, newsizes)
        cmd += " {0} /dev/null".format(netfile)
        sh(cmd)

    # Create liftOver chain file
    liftoverfile = pf + ".liftover.chain"
    if need_update((netfile, sortedchain), liftoverfile):
        cmd = "netChainSubset {0} {1} {2}".\
                format(netfile, sortedchain, liftoverfile)
        sh(cmd)

def print_chain(cid, tName, qName, qStrand, tSize, qSize, locs):
    chain = "chain"
    score = 1000
    tStrand = "+"
    tStart = min(x[0] for x in locs) - 1
    tEnd = max(x[1] for x in locs)
    qStart = min(x[2] for x in locs) - 1
    qEnd = max(x[3] for x in locs)
    headerline = " ".join(str(x) for x in (
         chain, score, tName, tSize, tStrand, tStart,
         tEnd, qName, qSize, qStrand, qStart, qEnd, cid
    ))
    print(headerline)
    for i in range(len(locs)):
        tb, te, qb, qe = locs[i]
        size, size2 = te - tb + 1, qe - qb + 1
        assert size == size2, "size not equal"
        if i == len(locs) - 1:
            print(size)
        else:
            dt = locs[i+1][0] - 1 - te
            if qStrand == "-":
                dq = qb - 1 - locs[i+1][3]
            else:
                dq = locs[i+1][2] - 1 - qe
            print("%d\t%d\t%d" % (size, dt, dq))
    print()

def bed2chain(args):
    from maize.formats.sizes import Sizes
    tdic = Sizes(args.tsize)
    qdic = Sizes(args.qsize)
    
    firstline = True
    cid0, tName0, qName0, srd0, locs = '', '', '', '', []
    for line in must_open(args.fi):
        line = line.rstrip("\n")
        if not line:
            continue
        tName, tStart, tEnd, srd, qName, qStart, qEnd, cid = line.split()
        tStart, tEnd, qStart, qEnd = int(tStart), int(tEnd), int(qStart), int(qEnd)
        if firstline:
            cid0, tName0, qName0, srd0 = cid, tName, qName, srd
            locs.append([tStart, tEnd, qStart, qEnd])
            firstline = False
        elif cid0 == cid:
            assert tName == tName0 and qName == qName0 and srd == srd0, "inconsistent info in chain"
            locs.append([tStart, tEnd, qStart, qEnd])
        else:
            print_chain(cid0, tName0, qName0, srd0, tdic.get_size(tName0), qdic.get_size(qName0), locs)
            cid0, tName0, qName0, srd0 = cid, tName, qName, srd
            locs = [[tStart, tEnd, qStart, qEnd]]
    print_chain(cid0, tName0, qName0, srd0, tdic.get_size(tName0), qdic.get_size(qName0), locs)
 
def chain2bed(args):
    chainFile = Chain(args.fi)
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
                print("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (qName, qb, qe, qSrd, tName, tb, te, cid))
            else:
                print("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s" % (tName, tb, te, qSrd, qName, qb, qe, cid))
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
    sp1.add_argument('--qry', action = 'store_true', help = 'use query coordinate system')
    sp1.set_defaults(func = chain2bed)
    
    sp1 = sp.add_parser("frombed", help = "convert from BED to chain file")
    sp1.add_argument('fi', help = 'input tsv file')
    sp1.add_argument('tsize', help = 'target size file')
    sp1.add_argument('qsize', help = 'query size file')
    sp1.set_defaults(func = bed2chain)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



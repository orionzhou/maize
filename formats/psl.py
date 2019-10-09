#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re
import logging

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import LineFile, must_open
from maize.formats.sizes import Sizes
from maize.utils.location import locAry2Str

class PslLine(object):

    def __init__(self, sline):
        args = sline.strip().split()
        self.nargs = len(args)
        self.matches = int(args[0])
        self.misMatches = int(args[1])
        self.repMatches = int(args[2])
        self.nCount = int(args[3])
        self.qNumInsert = int(args[4])
        self.qBaseInsert = int(args[5])
        self.tNumInsert = int(args[6])
        self.tBaseInsert = int(args[7])
        self.qstrand, self.gstrand = args[8], None
        #m = re.match(r"(?P<qs>[\+\-]?)(?P<gs>[\+\-])", self.qstrand)
        m = re.match(r"(?P<qs>[\+\-])(?P<gs>[\+\-])", self.qstrand)
        if m:
            self.qstrand, self.gstrand = m.group('qs'), m.group('gs')
        self.qName = args[9]
        self.qSize = int(args[10])
        self.qStart = int(args[11])
        self.qEnd = int(args[12])
##        if self.qstrand == "-":
##            self.qStart, self.qEnd = self.qSize - self.qEnd, \
##                    self.qSize - self.qStart
        self.tName = args[13]
        self.tSize = int(args[14])
        self.tStart = int(args[15])
        self.tEnd = int(args[16])
        self.blockCount = int(args[17])
        self.blockSizes = [int(x) for x in args[18].strip().split(',')[:-1]]
        self.qStarts = [int(x) for x in args[19].strip().split(',')[:-1]]
        self.tStarts = [int(x) for x in args[20].strip().split(',')[:-1]]
##        self.tStarts = [self.tSize - int(x) if self.strand == "-" \
##                else int(x) for x in args[20].strip().split(',')[:-1]]

    def __str__(self):
        args = [self.matches, self.misMatches, self.repMatches, \
                self.nCount, self.qNumInsert, self.qBaseInsert, \
                self.tNumInsert, self.tBaseInsert, self.qstrand, \
                self.qName, self.qSize, self.qStart, self.qEnd, \
                self.tName, self.tSize, self.tStart, self.tEnd, \
                self.blockCount, \
                ",".join(str(x) for x in self.blockSizes)+",", \
                ",".join(str(x) for x in self.qStarts)+",", \
                ",".join(str(x) for x in self.tStarts)+","]

        s = "\t".join(str(x) for x in args)
        return s

    def __getitem__(self, key):
        return getattr(self, key)

    @property
    def qspan(self):
        return self.qEnd - self.qStart

    @property
    def tspan(self):
        return self.tEnd - self.tStart

    @property
    def score(self):
        sizeMult = self._sizeMult

        return sizeMult * (self.matches + (self.repMatches >> 1)) - \
               sizeMult * self.misMatches - self.qNumInsert - self.tNumInsert

    @property
    def coverage(self):
        return 100 * (self.matches + self.misMatches + \
                self.repMatches + self.nCount) / self.qSize

    @property
    def swap(self):
        self.qName, self.qSize, self.tName, self.tSize = \
                self.tName, self.tSize, self.qName, self.qSize

        self.qStart, self.qEnd, self.tStart, self.tEnd = \
                self.tStart, self.tEnd, self.qStart, self.qEnd

        self.qStarts, self.tStarts = self.tStarts, self.qStarts

    @property
    def _sizeMult(self):
        """
        decide the size multiplier based on sequence space (protein/nucleotide)
        """
        return 3 if self._isProtein else 1

    @property
    def _isProtein(self):
        """
        check if blockSizes and scores are in the protein space or not
        """
        last = self.blockCount - 1
        return ((self.tEnd == self.tStarts[last] + 3 * self.blockSizes[last]) \
                and self.strand == "+") or \
                ((self.tStart == self.tSize - (self.tStarts[last] + 3 * self.blockSizes[last])\
                and self.strand == "-"))

    def _milliBad(self, ismRNA=False):
        """
        calculate badness in parts per thousand
        i.e. number of non-identical matches
        """
        sizeMult = self._sizeMult

        qAlnSize, tAlnSize = self.qspan * sizeMult, self.tspan
        alnSize = min(qAlnSize, tAlnSize)
        if alnSize <= 0:
            return 0

        sizeDiff = qAlnSize - tAlnSize
        if sizeDiff < 0:
            sizeDiff = 0 if ismRNA else -sizeDiff

        insertFactor = self.qNumInsert
        if not ismRNA:
            insertFactor += self.tNumInsert

        total = (self.matches + self.repMatches + self.misMatches) * sizeMult

        return (1000 * (self.misMatches * sizeMult + insertFactor + \
                round(3 * math.log(1 + sizeDiff)))) / total if total != 0 else 0

    def pct_id(self, simple=None):
        return 100.00 - self._milliBad(ismRNA=True) * 0.1 if not simple \
                else 100.00 * self.matches / (self.matches + self.misMatches)
                #else 100.00 * self.score / self.qSize

    def gffline(self, source="GMAP", type="match_part", primary_tag="Parent", \
           alt_score=None, suffix=".match", count=0):

        score = "." if type == "match_part" else "{0:.2f}".format(self.score)

        target = " ".join(str(x) for x in [self.qName, self.qStart, self.qEnd])

        attributes = [primary_tag + "=" + self.qName + suffix + str(count), "Target=" + target]
        if primary_tag == "ID":
            attributes.extend(["identity={0:.2f}".format(self.pct_id(simple=alt_score)),\
                    "coverage={0:.2f}".format(self.coverage)])
        attrs = ";".join(str(x) for x in attributes)

        line = "\t".join(str(x) for x in [self.tName, source, type, self.tStart, \
                self.tEnd, score, self.strand, ".", attrs])
        return line

    @property
    def bed12line(self):
        color = "255,0,0"
        self.blockStarts = ",".join([str(x - self.tStart) for x in self.tStarts])
        line = "\t".join(str(x) for x in (self.tName, self.tStart, self.tEnd, \
                self.qName, "{0:.2f}".format(self.pct_id()), self.strand, \
                self.tStart, self.tEnd, color, self.blockCount, \
                ",".join(str(bs) for bs in self.blockSizes), \
                self.blockStarts))
        return line

class Psl(LineFile):

    def __init__(self, filename=None):
        super(Psl, self).__init__(filename)

        import re

        self.mCounts = {}   # dict to hold match counts
        if not filename:
            return

        for line in must_open(filename):
            if not re.match(r'\d+', line[0]):
                continue
            self.append(PslLine(line))

    def trackMatches(self, id):
        self.mCounts[id] = self.mCounts.get(id, 0) + 1

    def getMatchCount(self, id):
        return self.mCounts[id]

def coordQ(args):
    sizes = Sizes(args.fs)
    for line in must_open(args.fi):
        if not re.match(r'\d+', line[0]):
            continue
        p = PslLine(line)
        qnames = p.qName.split("-")
        if len(qnames) == 3:
            x = p.qName
            p.qName, qosStart, qosEnd = qnames[0], int(qnames[1]), int(qnames[2])
            assert qosEnd-qosStart+1 == p.qSize
            cSize = sizes.get_size(p.qName)
            p.qStart += qosStart - 1
            p.qEnd += qosStart - 1
            if p.qstrand == "-":
                p.qStarts = [x + cSize - qosEnd for x in p.qStarts]
            else:
                p.qStarts = [x + qosStart - 1 for x in p.qStarts]
            p.qSize = cSize
        print(str(p))

def coordT(args):
    sizes = Sizes(args.fs)
    for line in must_open(args.fi):
        if not re.match(r'\d+', line[0]):
            continue
        p = PslLine(line)
        tnames = p.tName.split("-")
        if len(tnames) == 3:
            x = p.tName
            p.tName, tosStart, tosEnd = tnames[0], int(tnames[1]), int(tnames[2])
            assert tosEnd-tosStart+1 == p.tSize
            cSize = sizes.get_size(p.tName)
            p.tStart += tosStart - 1
            p.tEnd += tosStart - 1
            p.tStarts = [x + tosStart - 1 for x in p.tStarts]
            p.tSize = cSize
        print(str(p))

def psl2tsv(args):
    sMatch, sMisMatch, sGapOpen, sGapExtend = 2, -3, -5, -2
    print("\t".join('''qName qStart qEnd qSize strand
    tName tStart tEnd tSize
    alnLen match misMatch baseN qNumIns tNumIns qBaseIns tBaseIns ident score
    qLoc tLoc'''.split()))
    for line in must_open(args.fi):
        if not re.match(r'\d+', line[0]):
            continue
        p = PslLine(line)
        qName, qStart, qEnd, qSize, strand = p.qName, p.qStart, p.qEnd, p.qSize, p.qstrand
        tName, tStart, tEnd, tSize = p.tName, p.tStart, p.tEnd, p.tSize
        match, misMatch, baseN, qNumIns, tNumIns, qBaseIns, tBaseIns = \
                p.matches, p.misMatches, p.nCount, p.qNumInsert, p.tNumInsert, \
                p.qBaseInsert, p.tBaseInsert

        assert p.blockCount==len(p.tStarts), "unequal pieces"
        assert p.blockCount==len(p.qStarts), "unequal pieces"
        assert p.blockCount==len(p.blockSizes), "unequal pieces"
        match += p.repMatches
        alnLen = match + misMatch + baseN
        assert alnLen==sum(p.blockSizes), "block size error: %s %d %d" % (qId, alnLen, sum(blockSizes))
        assert alnLen+qBaseIns==qEnd-qStart, "%s: qLen error" % qId
        assert alnLen+tBaseIns==tEnd-tStart, "%s: tLen error" % qId

        qLocs, tLocs = [], []
        for i in range(p.blockCount):
            rtb, rte = p.tStarts[i]-tStart, p.tStarts[i]-tStart+p.blockSizes[i]
            rqb, rqe = 0, 0
            if strand == '-':
                rqb, rqe = p.qStarts[i]-(qSize-qEnd), p.qStarts[i]-(qSize-qEnd)+p.blockSizes[i]
            else:
                rqb, rqe = p.qStarts[i]-qStart, p.qStarts[i]+p.blockSizes[i]-qStart
            qLocs.append([rqb+1,rqe])
            tLocs.append([rtb+1,rte])

        score = match * sMatch + misMatch * sMisMatch
        numIns = qNumIns + tNumIns
        if numIns >= 1:
            score += sGapOpen + (numIns - 1) * sGapExtend
        ident = "%.03f" % (float(match)/(match+misMatch))
        print("\t".join(str(x) for x in [qName, qStart+1, qEnd, qSize, strand,
            tName, tStart+1, tEnd, tSize, alnLen, match, misMatch, baseN,
            qNumIns, tNumIns, qBaseIns, tBaseIns, ident, score,
            locAry2Str(qLocs), locAry2Str(tLocs)]))

def psl2bed(args):
    for line in must_open(args.fi):
        if not re.match(r'\d+', line[0]):
            continue
        p = PslLine(line)
        for i in range(p.blockCount):
            tbeg, tend = p.tStarts[i], p.tStarts[i] + p.blockSizes[i]
            qbeg, qend = 0, 0
            if p.qstrand == '-':
                qbeg, qend = p.qSize - p.qStarts[i] - p.blockSizes[i], p.qSize - p.qStarts[i]
            else:
                qbeg, qend = p.qStarts[i], p.qStarts[i] + p.blockSizes[i]
            tstr = "%s:%d-%d" % (p.tName, tbeg, tend)
            qstr = "%s:%d-%d" % (p.qName, qbeg, qend)
            if args.qry:
                print("%s\t%d\t%d\t%s\t%s" % 
                    (p.qName, qbeg, qend, p.qstrand, tstr))
            else:
                print("%s\t%d\t%d\t%s\t%s_%d_%d_%s" %
                    (p.tName, tbeg, tend, p.qstrand, qstr))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'psl utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("coordQ", help = "recover query chrom coordinates")
    sp1.add_argument('fi', help = 'input PSL')
    sp1.add_argument('fs', help = 'query genome size file (.sizes)')
    sp1.set_defaults(func = coordQ)
 
    sp1 = sp.add_parser("coordT", help = "recover target chrom coordinates")
    sp1.add_argument('fi', help = 'input PSL')
    sp1.add_argument('fs', help = 'target genome size file (.sizes)')
    sp1.set_defaults(func = coordT)
 
    sp1 = sp.add_parser("2bed", help = "convert to BED file")
    sp1.add_argument('fi', help = 'input PSL')
    sp1.add_argument('--qry', action = 'store_true', help = 'use query coordinate system')
    sp1.set_defaults(func = psl2bed)
    
    sp1 = sp.add_parser("2tsv", help = "convert to tsv file")
    sp1.add_argument('fi', help = 'input PSL')
    sp1.set_defaults(func = psl2tsv)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import numpy as np
import math
import pysam
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'unknown'
    )
    parser.add_argument(
        'org', nargs = "?", default = "HM034", help = 'genome to valid'
    )
    args = parser.parse_args()

    org = args.org
    dirw = op.join(os.environ['misc3'], "comp.sv.valid")
    fi = "%s/01_sv/%s.tbl" % (dirw, org)
    fo = "%s/02_pacbio/%s.tbl" % (dirw, org)
    
    fb = "%s/pacbio/%s_%s/15.bam" % (os.environ['misc3'], org, org)
    bam = pysam.Samfile(fb, "rb")
    
    fhi = open(fi, "r")
    fho = open(fo, "w")
    print >>fho, "\t".join(["id", "n1", "n2"])
    
    offset = 200
    next(fhi)
    for line in fhi:
        line = line.strip("\n")
        if line == "":
            break
        (id, tchr, tbeg, tend, tlen, tinfo, srd, qchr, qbeg, qend, qlen, qinfo) = line.split("\t")
        (qbeg, qend) = (int(qbeg), int(qend))
        iter = bam.fetch(qchr, qbeg, qend)
        (n1, n2) = (0, 0)
        for x in iter:
            if(x.pos+1 < qbeg - offset and x.aend > qbeg + offset):
                n1 += 1
            if(x.pos+1 < qend - offset and x.aend > qend + offset):
                n2 += 1
        print >>fho, "%s\t%d\t%d" % (id, n1, n2)
    bam.close()

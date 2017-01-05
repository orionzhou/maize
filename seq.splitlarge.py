#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import numpy as np
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'split large fasta record to even pieces'
    )
    parser.add_argument(
        'fi', help = 'input file (fasta)'
    )
    parser.add_argument(
        'fo', help = 'output file (fasta)'
    )
    parser.add_argument(
        '--opt', dest='opt', type=int, default=1, help='1: 100kb chunks (default); 2: 240 pieces'
    )
    args = parser.parse_args()

    fho = open(args.fo, "w")

    ary = []
    fhi = open(args.fi, "r")
    for seq in SeqIO.parse(fhi, "fasta") :
        ary.append(len(seq.seq))
    fhi.close()
    totalsize = sum(ary)
    
    if args.opt == 1:
        piecesize = 100000
    else:
        npieces = 10 * 24 
        piecesize = int((totalsize / npieces) / 100000) * 100000
    print("  total size: %d, size per piece: %d" % (totalsize, piecesize))
    
    fhi = open(args.fi, "r")
    for seq in SeqIO.parse(fhi, "fasta") :
        size = len(seq.seq)
        if(float(size) / piecesize > 1.3) :
            print("    splitting %s: %d" %(seq.id, size))
            ary = seq.id.split("-")
            [id, bbeg] = [ary[0], int(ary[1])]

            seqstr = str(seq.seq)
            nf = int(math.ceil(float(size) / piecesize))
            rcds = []
            for i in range(0, nf) :
                rbeg = i * piecesize
                rend = min((i+1) * piecesize, size)
                sseqstr = seqstr[rbeg:rend]
                sid = "%s-%d-%d" % (id, bbeg+rbeg, bbeg+rend-1)
                rcd = SeqRecord(Seq(sseqstr), id = sid, description = '')
                rcds.append(rcd)
                #print "      " + sid
            SeqIO.write(rcds, fho, "fasta")
        else:
            SeqIO.write(seq, fho, "fasta")
    fhi.close()
    fho.close()

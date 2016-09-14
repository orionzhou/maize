#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path
import argparse
from time import clock, time
from Seq import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='JJJ')
    parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, help='input')
    parser.add_argument('output', type=argparse.FileType('w'), default=sys.stdout, help='output')
    parser.parse_args()
    seq_db = "/home/youngn/zhoup/Data/misc3/spada.crp.Athaliana/01_genome/01_refseq.fa"
    seq_id = "Chr2"
    itv = 100000
    for i in range(1, 10):
        seq_beg = itv * i + 1
        seq_end = itv * i + 10
        c0, t0 = clock(), time()
        print seqret(seq_db, seq_id, seq_beg, seq_end, -1)
        print "process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0)
        

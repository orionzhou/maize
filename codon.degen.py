#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path
import argparse
from time import clock, time
from astropy.table import Table

if __name__ == "__main__":
#    parser = argparse.ArgumentParser(description='JJJ', prog="seqret")
#    parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, help='input')
#    parser.add_argument('output', type=argparse.FileType('w'), default=sys.stdout, help='output')
#    parser.parse_args()
    dirw = os.path.join(os.environ['misc2'], 'codon.degeneracy')
    fcod = os.path.join(dirw, 'codon.degeneracy.tbl')
    t = Table.read(fcod, format="ascii.tab")
    cdt = dict()
    for i in range(0, len(t)):
        cdt[t['codon'][i]] = [t['First'][i], t['Second'][i], t['Third'][i]]
#        print seqret(seq_db, seq_id, seq_beg, seq_end, -1)
#        print "process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0)
        

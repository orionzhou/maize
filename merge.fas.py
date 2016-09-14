#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import argparse

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
def read_cfg(fc):
    (orgs, fis) = ([], [])
    fhc = open(fc, "r")
    for line in fhc:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        (org, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            print os.access(fi, os.F_OK)
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
    fhc.close()
    return (orgs, fis)
def merge_seqs(fis, fids, fo):
    print "  merging input files to %s" % fo
    seqs = []
    for i in range(0,len(fids)):
        handle = 0
        if (fis[i].endswith(".gz")):
            handle = gzip.open(fis[i], "rb")
        else:
            handle = open(fis[i], "rU")
        seq_it = SeqIO.parse(handle, "fasta")
        handle.close

        seqs1 = [SeqRecord(rcd.seq, id = fids[i] + "|" + rcd.id,
            description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'merge multiple fasta files and update IDs'
    )
 #   parser.add_argument('program', type=str, help='progname')
    parser.add_argument(
        'cfgfile', help = 'config file (a text file with identifier followed by the absolute path of fasta in each line)'
    )
    parser.add_argument(
        'outfile', help = 'output fasta'
    )
    args = parser.parse_args()

    (fc, fo) = (args.cfgfile, args.outfile)
    (fids, fis) = read_cfg(fc)
    
    print "%d input fasta detected" % len(fis)
    print "output file: %s" % fo
    merge_seqs(fis, fids, fo)

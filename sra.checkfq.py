#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'check fastq files and generate 11.tsv'
    )
    parser.add_argument(
        'dirw', help = 'working directory'
    )
    parser.add_argument(
        'paired', nargs = "?", default=False, help = 'paired-end (default: no)'
    )
    args = parser.parse_args()
   
    dirw, paired = args.dirw, args.paired
    os.chdir(dirw)

    assert op.isfile("01.srr.tsv"), "no 01.srr.tsv in %s" % dirw
    ary = np.genfromtxt("01.srr.tsv", names = True, dtype = None, delimiter = "\t")
    cols = list(ary.dtype.names)
    d05 = op.abspath("05.reads")
    assert op.isdir(d05), "no 05.reads in %s" % dirw
    fho = open("06.tsv", "w")
    if paired:
        print >>fho, "\t".join([cols[0], "dirf", "read1", "read2"]+cols[1:len(cols)])
    else:
        print >>fho, "\t".join([cols[0], "dirf", "read1"]+cols[1:len(cols)])
    for row in ary:
        row = list(row)
        rid = row[0]
        read1, read2 = "%s_1.fastq.gz" % rid, "%s_2.fastq.gz" % rid
        f1 = "%s/%s" % (d05, read1)
        f2 = "%s/%s" % (d05, read2)
        if paired:
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            print >>fho, "\t".join([rid, d05, read1, read2] + row[1:len(row)])
        else:
            assert op.isfile(f1), "%s not there" % f1
            print >>fho, "\t".join([rid, d05, read1] + row[1:len(row)])
    fho.close()

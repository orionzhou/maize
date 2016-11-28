#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'check fastq files and generate 00.1.read.tsv'
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

    fi, fo = "00.0.srr.tsv", "00.1.read.tsv"
    assert op.isfile(fi), "no %s in %s" % (fi, dirw)
    ary = np.genfromtxt(fi, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    d05 = op.abspath("05.reads")
    assert op.isdir(d05), "no 05.reads in %s" % dirw
    fho = open(fo, "w")
    if paired:
        print >>fho, "\t".join(cols + ["ReadFile1", "ReadFile2"])
    else:
        print >>fho, "\t".join(cols + ["ReadFile"])
    for row in ary:
        row = list(row)
        sid = row[0]
        read1, read2 = "%s_1.fastq.gz" % sid, "%s_2.fastq.gz" % sid
        f1 = "%s/%s" % (d05, read1)
        f2 = "%s/%s" % (d05, read2)
        if paired:
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            print >>fho, "\t".join(row + [f1, f2])
        else:
            assert op.isfile(f1), "%s not there" % f1
            print >>fho, "\t".join(row + [f1])
    fho.close()

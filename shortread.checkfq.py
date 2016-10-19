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
        'pid', help = 'PID (under $data/misc2/rnaseq)'
    )
    args = parser.parse_args()
   
    pid = args.pid
    dirw = op.join(os.environ['misc2'], 'rnaseq', pid)
    os.chdir(dirw)

    assert op.isfile("01.srr.tsv"), "no 01.srr.tsv in %s" % dirw
    ary = np.genfromtxt("01.srr.tsv", names = True, dtype = None, delimiter = "\t")
    cols = list(ary.dtype.names)
    d05 = op.abspath("05.reads")
    assert op.isdir(d05), "no 05.reads in %s" % dirw
    fho = open("11.tsv", "w")
    print >>fho, "\t".join([cols[0], "PE", "dirf", "read1", "read2"]+cols[1:len(cols)])
    for row in ary:
        row = list(row)
        rid = row[0]
        read1, read2 = "%s_1.fastq.gz" % rid, "%s_2.fastq.gz" % rid
        f1 = "%s/%s" % (d05, read1)
        f2 = "%s/%s" % (d05, read2)
        assert op.isfile(f1), "%s not there" % f1
        pe = '0'
        if op.isfile(f2):
            pe = '1'
        else:
            read2 = ""
        print >>fho, "\t".join([rid, pe, d05, read1, read2] + row[1:len(row)])
    fho.close()

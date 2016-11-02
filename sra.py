#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import numpy as np
import sys
import time
import argparse
import itertools

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'download SRA files and extract fastq sequences'
    )
    parser.add_argument(
        'dirw', default = "/scratch.global/zhoux379/shortread/test", help = 'working directory (default: /scratch.global/zhoux379/shortread/test)'
    )
    parser.add_argument(
            'sralist', default = "/scratch.global/zhoux379/shortread/test/01.srr.tsv", help = 'list file of SRA IDs (default: /scratch.global/zhoux379/shortread/test/01.srr.tsv)'
    )
    args = parser.parse_args()

    dirw, sralist = args.dirw, args.sralist
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(sralist), "%s is not a file" % sralist 
    ary = np.genfromtxt(sralist, names = True, dtype = None, delimiter = "\t")
    cols = list(ary.dtype.names)
    if not op.isdir("03.sra"): os.makedirs("03.sra")
    if not op.isdir("05.reads"): os.makedirs("05.reads")
    d03 = "03.sra" #op.abspath("03.sra")
    d05 = "05.reads" #op.abspath("05.reads")
    fhc = open("04.fastqdump.sh", "w")
    os.chdir("03.sra")
    fnames = []
    for row in ary:
        rid = row[0]
        print rid
        cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra" % (rid[0:3], rid[0:6], rid, rid)
        #os.system(cmd)
        f_sra = "%s.sra" % rid
        assert op.isfile(f_sra), "%s is not there" % f_sra
        fnames.append(f_sra)
        print >>fhc, "fastq-dump --gzip --split-files -outdir %s %s/%s" % (d05, d03, f_sra)
    fname_str = " ".join(fnames)
    fhc.close()

    name = 'sra'
    os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
    assert op.isfile(name), "no %s in pbs" % name
    os.system("qsub %s -v DIR=%s" % (name, dirw))

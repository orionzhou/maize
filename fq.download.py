#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import numpy as np
from string import Template
import sys
import time
import argparse
import configparser
from colorama import init, Fore, Back, Style

def fq_download(dirw, ilist, olist, diro, paired):
    assert op.isfile(ilist), "%s is not a file" % ilist 
    if not op.isdir(diro): os.makedirs(diro)
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    if ary.size <= 1: ary = [ary.tolist()]
    
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + ["ReadFile1", "ReadFile2"])+"\n")
    else:
        fho.write("\t".join(cols + ["ReadFile"])+"\n")
    init()
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        print(Fore.GREEN + "working on %s" % sid + Style.RESET_ALL)
        cmd = "%s --gzip --split-files -outdir %s %s\n" % \
                (fastq_dump, diro, sid)
        os.system(cmd)
        read1, read2 = "%s_1.fastq.gz" % sid, "%s_2.fastq.gz" % sid
        f1 = "%s/%s" % (diro, read1)
        f2 = "%s/%s" % (diro, read2)
        if paired:
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            fho.write("\t".join(row + [f1, f2])+"\n")
        else:
            assert op.isfile(f1), "%s not there" % f1
            fho.write("\t".join(row + [f1])+"\n")
    fho.close()
def fq_download_check(dirw, ilist, olist, diro, paired):
    assert op.isfile(ilist), "%s is not a file" % ilist 
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    if ary.size <= 1: ary = [ary.tolist()]
    
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + ["ReadFile1", "ReadFile2"])+"\n")
    else:
        fho.write("\t".join(cols + ["ReadFile"])+"\n")
    init()
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        read1, read2 = "%s_1.fastq.gz" % sid, "%s_2.fastq.gz" % sid
        f1 = "%s/%s" % (diro, read1)
        f2 = "%s/%s" % (diro, read2)
        if paired:
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            fho.write("\t".join(row + [f1, f2])+"\n")
        else:
            assert op.isfile(f1), "%s not there" % f1
            fho.write("\t".join(row + [f1])+"\n")
    fho.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'Download and extract Fastq sequences from NCBI-SRA'
    )
    parser.add_argument(
            'config', nargs = '?', default = "config.ini", \
                    help = 'config file (default: config.ini)'
    )
    parser.add_argument(
            '--check', action = "store_true", \
                    help = 'run the script in check mode (default: no)'
    )
    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    cfg = cfg['fastq_download']
    dirw, ilist, olist = \
            cfg['dirw'], cfg['ilist'], cfg['olist']
    diro = cfg['outdir']
    paired = cfg.getboolean('paired')
    fastq_dump = cfg['fastq_dump']
    
    if args.check:
        fq_download_check(dirw, ilist, olist, diro, paired)
        sys.exit(0)
    fq_download(dirw, ilist, olist, diro, paired)


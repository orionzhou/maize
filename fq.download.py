#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column

from maize.apps.base import AttrDict, str2bool, eprint, sh, mkdir, which
from maize.formats.base import must_open

def fq_download(ilist, diro, paired):
    assert op.isfile(ilist), "%s is not a file" % ilist 
    if not op.isdir(diro): os.makedirs(diro)
    t = Table.read(ilist, format = 'ascii.tab')
    for i in range(len(t)):
        sid = t['sid'][i]
        flag = False
        if paired:
            f1 = "%s/%s_1.fastq.gz" % (diro, sid)
            f2 = "%s/%s_2.fastq.gz" % (diro, sid)
            if op.isfile(f1) and op.isfile(f2):
                flag = True
        else:
            f1 = "%s/%s.fastq.gz" % (diro, sid)
            if op.isfile(f1):
                flag = True
        splitflag = '--split-files' if paired else ''
        cmd = "%s --gzip %s -outdir %s %s" % \
                ('fastq-dump', splitflag, diro, sid)
        if flag:
            logging.debug("%s already done - skip" % sid)
        else:
            logging.debug("working on %s" % sid)
            sh(cmd)
def fq_download_check(ilist, diro, paired):
    assert op.isfile(ilist), "%s is not a file" % ilist 
    if not op.isdir(diro): os.makedirs(diro)
    t = Table.read(ilist, format = 'ascii.tab')
    for i in range(len(t)):
        sid = t['sid'][i]
        if paired:
            read1, read2 = "%s_1.fastq.gz" % sid, "%s_2.fastq.gz" % sid
            f1 = "%s/%s" % (diro, read1)
            f2 = "%s/%s" % (diro, read2)
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
        else:
            f1 = "%s/%s.fastq.gz" % (diro, sid)
            assert op.isfile(f1), "%s not there" % f1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Download and extract Fastq sequences from NCBI-SRA'
    )
    parser.add_argument(
            'ilist', help = 'list of SRA numbers'
    )
    parser.add_argument(
            'outdir', help = 'output directory'
    )
    parser.add_argument(
            '--paired', action = 'store_true', 
            help = 'paired-end reads'
    )
    parser.add_argument(
            '--check', action = "store_true",
             help = 'run script in check mode'
    )
    args = parser.parse_args()
    ilist, diro, paired = args.ilist, args.outdir, args.paired
    if args.check:
        fq_download_check(ilist, diro, paired)
        sys.exit(0)
    fq_download(ilist, diro, paired)


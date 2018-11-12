#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column

from maize.formats.base import must_open
from maize.apps.base import AttrDict, str2bool, eprint, sh, mkdir, which


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'rename fastq files'
    )
    parser.add_argument('diri', help = 'FROM directory path')
    parser.add_argument('diro', help = 'TO directory path')
    parser.add_argument('list', help = 'read list (sid Readfile1[ Readfile2])')
    args = parser.parse_args()
    
    #dirw = '/home/springer/zhoux379/projects/3rnaseq/cache'
    #diri = '/scratch.global/tenders/3RNA_0418/fastq_gz_files'
    #diro = '10.fastq'
    diri, diro, ilist = args.diri, args.diro, args.list
    
    if not op.isfile(ilist):
        logging.error("cannot read %s" % ilist)
        sys.exit()
    if not op.isdir(diro):
        mkdir(diro)
    
    t = Table.read(ilist, format = 'ascii.tab')
    paired = False
    if 'Readfile2' in t.colnames:
        paired = True
    tag = "paired" if paired else "single"
    logging.debug("proceed as %s-end reads" % tag)

    for i in range(len(t)):
        if paired:
            sid, fq1, fq2 = t['sid'][i], t['Readfile1'][i], t['Readfile2'][i]
            fq1, fq2 = op.join(diri, fq1), op.join(diri, fq2)
            assert op.isfile(fq1), "cannot access fq1: %s" % fq1
            assert op.isfile(fq2), "cannot access fq2: %s" % fq2
            cmd = "cp %s %s/%s_1.fq.gz" % (fq1, diro, sid)
            sh(cmd)
            cmd = "cp %s %s/%s_2.fq.gz" % (fq2, diro, sid)
            sh(cmd)
        else:
            sid, fq = t['sid'][i], t['Readfile'][i]
            fq = op.join(diri, fq)
            assert op.isfile(fq), "cannot access fq: %s" % fq
            cmd = "cp %s %s/%s.fq.gz" % (fq, diro, sid)
            sh(cmd)


#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import argparse
import configparser
from string import Template
from colorama import init, Fore, Back, Style
import pysam

def check_bam(fbam):
    exist_bam = 1
    try:
        bam = pysam.AlignmentFile(fbam, "rb")
    except:
        exist_bam = 0
    return exist_bam
def check_sam(fsam):
    exist_sam = 1
    try:
        sam = pysam.AlignmentFile(fsam, "r")
    except:
        exist_sam = 0
    return exist_sam

def bwa(dirw, ilist, olist, do1, paired, db_bwa,
        bwa, snpEff, samtools, parallel,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1 = "41.1.bwa.sh"
    fho1 = open(fo1, "w")
    for diro in [do1]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        fsam = "%s/%s.sam" % (do1, sid)
        exist_fsam = check_sam(fsam)
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            fho1.write("%s mem -t %d -M %s %s %s > %s/%s.PE.sam\n" % \
                    (bwa, pbs_ppn, db_bwa, f1p, f2p, do1, sid))
            fho1.write("%s mem -t %d -P -S -M %s %s %s > %s/%s.SE.sam\n" % \
                    (bwa, pbs_ppn, db_bwa, f1u, f2u, do1, sid))
            fho.write("java -jar picard.jar MergeSamFiles I=%s/%s.PE.sam \
                    I=%s/%s.SE.sam O=%s/%s.sam\n" % \
                    (do1, sid, do1, sid, do1, sid))
        else:
            fr, rc, ft, rrc = row[5:9]
            fho1.write("%s mem -t %d -P -S -M %s %s > %s.sam\n" % \
                    (bwa, pbs_ppn, db_bwa, ft, sid))
    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("module load picard/2.3.0")
    cmds.append("module load gatk/3.7.0")
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo1))
    cmd = "\n".join(cmds)
    
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": cmd
    }
    fo = "41.pbs"
    fho = open(fo, "w")
    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    fho.write(src.substitute(temdict))
    
    init()
    print(Fore.GREEN)
    print("A job script has been generated: %s" % fo)
    print("Please check, make necessary changes , then type:")
    print(Fore.RED + "qsub %s" % fo)
    print(Style.RESET_ALL)
def bwa_check(dirw, ilist, olist, diro):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    fho.write("\t".join(cols + ["HtseqFile"])+"\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        fhts = "%s/%s.txt" % (diro, sid)
        assert op.isfile(fhts), "%s not there" % fhts
        fho.write("\t".join(row + [fhts]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Quantify gene expression using htseq'
    )
    parser.add_argument(
            'config', nargs = '?', default = "config.ini", help = 'config file (default: config.ini)'
    )
    parser.add_argument(
            '--check', action = "store_true", help = 'run the script in check mode (default: no)'
    )
    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    cfg = cfg['bwa']
    dirw, ilist, olist, do1 = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['outdir1']
    paired = cfg.getboolean('paired')
    annotation = cfg['annotation']
    db_bwa = cfg['db_bwa']
    bwa, snpEff, samtools, parallel = \
            cfg['bwa'], cfg['snpEff'], cfg['samtools'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if args.check:
        bwa_check(dirw, ilist, olist, diro)
        sys.exit(0)
    bwa(dirw, ilist, olist, do1, paired, db_bwa,
            bwa, snpEff, samtools, parallel,
            pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email)


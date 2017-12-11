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

def run_htseq(dirw, ilist, olist, diro, paired, srd, annotation,
        htseq, samtools, parallel,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1 = "31.1.htseq.sh"
    fho1 = open(fo1, "w")
    for diro in [diro]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        fsen = "%s/%s.txt" % (diro, sid)
        fant = "%s/%s.as.txt" % (diro, sid)
        if paired:
            fbam = row[15]
            if not op.isfile(fsen) or os.stat(fsen).st_size == 0:
                fho1.write("%s view -f 1 -F 256 %s | %s -r pos -s %s \
                        -t exon -i gene_id -m union -a 20 - %s > %s/%s.txt\n" % 
                        (samtools, fbam, htseq, 'reverse', annotation, diro, sid))
            if not op.isfile(fant) or os.stat(fant).st_size == 0:
                fho1.write("%s view -f 1 -F 256 %s | %s -r pos -s %s \
                        -t exon -i gene_id -m union -a 20 - %s > %s/%s.as.txt\n" % 
                        (samtools, fbam, htseq, 'yes', annotation, diro, sid))
        else:
            fbam = row[9]
            if not op.isfile(fsen) or os.stat(fsen).st_size == 0:
                fho1.write("%s view -F 256 %s | %s -r pos -s no \
                        -t exon -i gene_id -m union -a 20 - %s > %s/%s.txt\n" % 
                        (samtools, fbam, htseq, annotation, diro, sid))

    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("module load python2")
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo1))
    cmd = "\n".join(cmds)
    
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": cmd
    }
    fo = "31.pbs"
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
def htseq_check(dirw, ilist, olist, diro):
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
            'config', nargs = '?', default = "config.ini", 
            help = 'config file (default: config.ini)'
    )
    parser.add_argument(
            '--check', action = "store_true", 
            help = 'run the script in check mode (default: no)'
    )
    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    cfg = cfg['htseq']
    dirw, ilist, olist, diro = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['outdir']
    paired = cfg.getboolean('paired')
    srd, annotation = cfg['stranded'], cfg['annotation']
    htseq, samtools, parallel = cfg['htseq'], cfg['samtools'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if args.check:
        htseq_check(dirw, ilist, olist, diro)
        sys.exit(0)
    run_htseq(dirw, ilist, olist, diro, paired, srd, annotation,
            htseq, samtools, parallel,
            pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email)


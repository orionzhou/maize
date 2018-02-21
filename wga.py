#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import pysam

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob

def prepare(cfg):
    cfg = cfg['prepare']
    dirw, qry, tgt = cfg['dirw'], cfg['qry'], cfg['tgt']
    qry_fas, tgt_fas = cfg['qry_fas'], cfg['tgt_fas']
    tmpdir = cfg['temp_dir']
    
    if not op.isdir(dirw):
        logging.debug("making directory: %s" % dirw)
        mkdir(dirw)
    os.chdir(dirw)
    subdirs = ['01_tgt_genome', '02_qry_genome']
    for subdir in subdirs:
        if not op.isdir(subdir):
            logging.debug("making directory: %s" % subdir)
            mkdir(subdir)

    if op.isfile(tgt_fas):
        fo = "raw.fa"
        if tgt_fas.endswith(".fa.gz") or tgt_fas.endswith(".fas.gz"):
            fo = "raw.fa.gz"
        sh("ln -sf %s 01_tgt_genome/%s" % (tgt_fas, fo))
    else:
        logging.error("%s not exist" % qry_fas)
    if op.isfile(qry_fas):
        fo = "raw.fa"
        if qry_fas.endswith(".fa.gz") or qry_fas.endswith(".fas.gz"):
            fo = "raw.fa.gz"
        sh("ln -sf %s 02_qry_genome/%s" % (qry_fas, fo))
    else:
        logging.error("%s not exist" % qry_fas)

    sh("genome %s/%s fasta --norename" % (dirw, "01_tgt_genome"))
    sh("genome %s/%s blat" % (dirw, "01_tgt_genome"))
    sh("genome %s/%s fasta --norename" % (dirw, "02_qry_genome"))

def run_pblat(cfg):
    cfg = cfg['blat']
    dirw, jobpre, diro = \
            cfg['dirw'], cfg['job_prefix'], cfg['outdir']
    qry, tgt = cfg['qry'], cfg['tgt']
    samtools, bcftools, parallel = \
            cfg['samtools'], cfg['bcftools'], cfg['parallel']

    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    if op.isdir(dirj):
        os.system("rm -rf %s" % dirj)
    for do in [diro, dirj]:
        if not op.isdir(do): 
            os.makedirs(do)
    fj = "%s.sh" % jobpre
    fhj = open(fj, "w")
    i = 1
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        gt = row[3]
        if paired:
            fbam = row[15]
        else:
            fbam = row[9]
        pre = "%s/%s" % (diro, sid)
        cmds = [
            "mkdir %s" % pre,
            "bam2bed.py %s %s.1.bed" % (fbam, pre),
            "sort -T %s -k1,1 -k2,2n %s.1.bed > %s.2.sorted.bed" % (pre, pre, pre),
            "intersectBed -wa -wb -a %s.2.sorted.bed -b %s > %s.3.bed" % (pre, target_vcf, pre),
            "sort -T %s -k4,4 -k1,1 -k2,2n %s.3.bed > %s.4.sorted.bed" % (pre, pre, pre),
            "bed.ase.py %s.4.sorted.bed %s.5.tsv %s.6.bed" % (pre, pre, pre),
            "sort -T %s -k1,1 -k2,2n %s.6.bed > %s.7.sorted.bed" % (pre, pre, pre),
            "intersectBed -wa -wb -a %s -b %s.7.sorted.bed > %s.8.bed" % (gene_bed, pre, pre),
            "bed.ase.sum.py %s.5.tsv %s.8.bed %s.tsv" (pre, pre, pre),
            "rm %s.[1-8].*" % pre,
            "rm -rf %s" % pre,
        ]
        fo = "%s/%03d.sh" % (dirj, i)
        fho = open(fo, "w")
        fho.write("\n".join(cmds) + "\n")
        fho.close()
        i += 1
        if not op.isfile("%s.bed" % pre) or not op.isfile("%s.tsv" % pre):
            fhj.write("bash %s\n" % fo)
    fhj.close()

    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    
    cmds = [
        "cd %s" % dirw,
        "%s -j %s < %s" % (parallel, pbs_ppn, fj)
    ]
    fjob = "%s.pbs" % jobpre
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": "\n".join(cmds)
    }
    fho = open(fjob, "w")
    fho.write(src.substitute(temdict))
    fho.close()

    init()
    print(Fore.GREEN)
    print("One job script has been generated: %s" % fjob)
    print("Please check, make necessary changes, then type:")
    print(Fore.RED + "qsub %s" % fjob)
    print(Style.RESET_ALL)

if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'whole genome alignment pipeline'
    )
    parser.add_argument('config', nargs = "?", default = "config.ini", help = 'config file')
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("prepare", help = "prepare sequences for align") 
    sp1.set_defaults(func = prepare)
    
    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    if args.command:
        args.func(cfg)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


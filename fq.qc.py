#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import argparse
import configparser
from fadapa import Fadapa
from string import Template
from colorama import init, Fore, Back, Style

def parse_fastqc(fqc):
    assert op.isfile(fqc), "%s not exist" % fqc
    qc = Fadapa(fqc)
    r = dict()
    for ary in qc.clean_data('Basic Statistics'):
        r[ary[0]] = ary[1]
    return r

def fq_trim(dirw, ilist, olist, do1, do2, do3, 
        paired, f_adp, fastqc, trimm, parallel,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(f_adp), "%s not exist" % f_adp 
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2, fo3 = "12.1.fastqc.sh", "12.2.trim.sh", "12.3.fastqc.sh"
    fho1, fho2, fho3 = open(fo1, "w"), open(fo2, "w"), open(fo3, "w")
    assert op.isfile(trimm), "%s is not there" % trimm
    for diro in [diro1, diro2, diro3]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        if paired:
            f1, f2 = row[5:7]
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            f11 = "%s/%s_1.PE.fastq.gz" % (do2, sid)
            f12 = "%s/%s_1.SE.fastq.gz" % (do2, sid)
            f21 = "%s/%s_2.PE.fastq.gz" % (do2, sid)
            f22 = "%s/%s_2.SE.fastq.gz" % (do2, sid)
            fho1.write("%s -o %s --extract -f fastq %s %s\n" % \
                    (fastqc, do1, f1, f2))
            fho2.write("java -Xmx2500M -jar %s PE -threads 4 \
                    %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10:8:no \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n" % \
                    (trimm, f1, f2, f11, f12, f21, f22, f_adp))
            fho3.write("%s -o %s --extract -f fastq %s %s %s %s\n" % \
                    (fastqc, do3, f11, f12, f21, f22))
        else:
            f1 = row[5]
            assert op.isfile(f1), "%s not there" % f1
            fho1.write("%s -o %s --noextract -f fastq %s\n" % \
                    (fastqc, diro1, f1))
            fo = "%s/%s.fastq.gz" % (do2, sid)
            fob = "%s/%s_1.fastq.gz" % (do2, sid)
            if op.isfile(fob):
                os.system("mv %s %s" % (fob, fo))
            fho2.write("java -Xmx2500M -jar %s SE -threads 4 \
                    %s %s ILLUMINACLIP:%s:2:30:10:8:no \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n" % \
                    (trimm, f1, fo, f_adp))
            fho3.write("%s -o %s --noextract -f fastq %s\n" % \
                    (fastqc, do3, fo))
    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo1))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo2))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo3))
    cmd = "\n".join(cmds)
    
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": cmd
    }
    fo = "12.pbs"
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
def fq_trim_check(dirw, ilist, olist, do1, do2, do3, paired):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + [ \
                "ReadPairCount", "TrimmedReadFile1Paired", \
                "TrimmedReadFile1Unpaired", "TrimmedReadFile2Paired", \
                "TrimmedReadFile2Unpaired", "RetainedReadPairCount", \
                "RetainedRead1Count", "RetainedRead2Count"]) + "\n")
    else:
        fho.write("\t".join(cols + [ \
                "ReadCount", "TrimmedReadFile", "TrimmedReadCount"]) + "\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        if paired:
            f1, f2 = row[5:7]
            f1p = "%s_1.PE.fastq.gz" % (sid)
            f1s = "%s_1.SE.fastq.gz" % (sid)
            f2p = "%s_2.PE.fastq.gz" % (sid)
            f2s = "%s_2.SE.fastq.gz" % (sid)
            p1p, p1s, p2p, p2s = [op.abspath("%s/%s" % (do2, x)) for x in [f1p, f1s, f2p, f2s]]
            assert op.isfile(p1p), "%s not exist" % p1p
            assert op.isfile(p1s), "%s not exist" % p1s
            assert op.isfile(p2p), "%s not exist" % p2p
            assert op.isfile(p2s), "%s not exist" % p2s
            pre1 = op.basename(f1).rstrip(".fastq.gz")
            pre2 = op.basename(f1).rstrip(".fastq.gz")
            q11 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre1)
            q12 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre2)
            r11, r12 = parse_fastqc(q11), parse_fastqc(q12)
            rc11, rc12 = r11['Total Sequences'], r12['Total Sequences']
            assert rc11 == rc12, "%s: read1 %d != read2 %d" % \
                    (sid, rc11, rc12)
            q21 = "%s/%s_1.PE_fastqc/fastqc_data.txt" % (do3, sid)
            q22 = "%s/%s_2.PE_fastqc/fastqc_data.txt" % (do3, sid)
            r21, r22 = parse_fastqc(q21), parse_fastqc(q22)
            rc21, rc22 = r21['Total Sequences'], r22['Total Sequences']
            assert rc21 == rc22, "%s: read1 %d != read2 %d" % \
                    (sid, rc21, rc22)
            q31 = "%s/%s_1.SE_fastqc/fastqc_data.txt" % (do3, sid)
            q32 = "%s/%s_2.SE_fastqc/fastqc_data.txt" % (do3, sid)
            r31, r32 = parse_fastqc(q31), parse_fastqc(q32)
            rc31, rc32 = r31['Total Sequences'], r32['Total Sequences']
            fho.write("\t".join(row + [ \
                    rc11, p1p, p1s, p2p, p2s, rc21, rc31, rc32]) + "\n")
        else:
            f1 = row[5]
            p1 = op.abspath("%s/%s.fastq.gz" % (do2, sid))
            assert op.isfile(p1), "%s not exist" % p1
            pre1 = op.basename(f1).rstrip(".fastq.gz")
            q1 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre1)
            q2 = "%s/%s_fastqc/fastqc_data.txt" % (do3, sid)
            r1, r2 = parse_fastqc(q1), parse_fastqc(q2)
            rc1, rc2 = r1['Total Sequences'], r2['Total Sequences']
            fho.write("\t".join(row + [rc1, p1, rc2]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Trimming and QC on a list of Fastq files'
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
    cfg = cfg['fastq_trim']
    dirw, ilist, olist, do1, do2, do3 = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], \
            cfg['outdir1'], cfg['outdir2'], cfg['outdir3']
    paired = cfg.getboolean('paired')
    f_adp, fastqc, trimm, parallel = \
            cfg['adapter'], cfg['fastqc'], cfg['trimmomatic'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if args.check:
        fq_trim_check(dirw, ilist, olist, do1, do2, do3, paired)
        sys.exit(0)
    fq_trim(dirw, ilist, olist, do1, do2, do3,
            paired, f_adp, fastqc, trimm, parallel,
            pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email)


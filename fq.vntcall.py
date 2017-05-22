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

def vntcall(dirw, ilist, olist, do1, do2, paired, ref_gatk,
        gatk, parallel, temp_dir,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2, fo3, fo4 = "50.1.rtc.sh", "50.2.ir.sh", "50.3.hc.sh", \
            "50.4.filter.sh"
    fho1, fho2, fho3, fho4 = [open(x, "w") for x in (fo1, fo2, fo3, fo4)]
    for diro in [do1, do2]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    bams = dict()
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        spec, tiss, geno = row[1:4]
        if geno not in bams:
            bams[geno] = []
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            bam = row[15]
        else:
            fr, rc, ft, rrc = row[5:9]
            bam = row[9]
        bams[geno].append("-I %s" % bam)
    jgatk = "java -jar %s" % gatk
    for geno in bams.keys(): 
        bamstr = " ".join(bams[geno])
        pre1 = "%s/%s" % (do1, geno)
        fho1.write("%s -T RealignerTargetCreator -nt %s \
                -R %s %s -o %s.list\n" % \
                (jgatk, pbs_ppn, ref_gatk, bamstr, pre1))
        fho2.write("%s -T IndelRealigner -R %s \
                %s -targetIntervals %s.list -o %s.bam\n" % \
                (jgatk, ref_gatk, bamstr, pre1, pre1))
        pre2 = "%s/%s" % (do2, geno)
        #fho3.write("%s -T UnifiedGenotyper -nt 8 -nct 3 -R %s \
        #        -I %s.bam -o %s.vcf\n" % \
        #        (jgatk, ref_gatk, pre1, pre2))
        fho3.write("%s -T HaplotypeCaller -nct %s -R %s \
                -I %s.bam -o %s.vcf\n" % \
                (jgatk, pbs_ppn, ref_gatk, pre1, pre2))
        fho4.write("%s -T SelectVariants -R %s \
                -V %s.vcf -selectType SNP -o %s.snp.raw.vcf\n" % \
                (jgatk, ref_gatk, pre2, pre2))
        fho4.write("%s -T SelectVariants -R %s \
                -V %s.vcf -selectType INDEL -o %s.indel.raw.vcf\n" % \
                (jgatk, ref_gatk, pre2, pre2))
        fho4.write("%s -T VariantFiltration --logging_level ERROR -R %s \
                -V %s.snp.raw.vcf --filterExpression 'QD < 2.0 \
                || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 \
                || ReadPosRankSum < -8.0 || SOR > 4.0' \
                --filterName 'basic_snp_filter' -o %s.snp.vcf\n" % \
                (jgatk, ref_gatk, pre2, pre2))
        fho4.write("%s -T VariantFiltration --logging_level ERROR -R %s \
                -V %s.indel.raw.vcf --filterExpression 'QD < 2.0 \
                || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
                --filterName 'basic_indel_filter' -o %s.indel.vcf\n" % \
                (jgatk, ref_gatk, pre2, pre2))
        fho4.write("%s -T CombineVariants \
                -R %s --assumeIdenticalSamples \
                -V %s.snp.vcf -V %s.indel.vcf \
                -o %s.filtered.vcf\n" \
                % (jgatk, ref_gatk, pre2, pre2, pre2))
        fho4.write("vcf2tsv.py %s.filtered.vcf %s.filtered.tsv\n" \
                % (pre2, pre2))
    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("module load java/jdk1.8.0_45")
    cmds.append("export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir)
    ##cmds.append("module load gatk/3.7.0")
    cmds.append("bash %s" % fo1)
    cmds.append("bash %s" % fo2)
    cmds.append("bash %s" % fo3)
    cmds.append("bash %s" % fo4)
    #cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo4))
    cmd = "\n".join(cmds)
    
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": cmd
    }
    fo = "50.pbs"
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
def vntcall_check(dirw, ilist, olist, do1, do2):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    fho.write("\t".join(cols + ["BAM"])+"\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        fhts = "%s/%s.txt" % (do1, sid)
        assert op.isfile(fhts), "%s not there" % fhts
        fho.write("\t".join(row + [fhts]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'Variant calling'
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
    cfg = cfg['vntcall']
    dirw, ilist, olist, do1, do2 = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['outdir1'], \
            cfg['outdir2']
    paired = cfg.getboolean('paired')
    annotation = cfg['annotation']
    ref_gatk = cfg['ref_gatk']
    temp_dir = cfg['temp_dir']
    parallel, gatk = cfg['parallel'], cfg['gatk']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if args.check:
        vntcall_check(dirw, ilist, olist, do1, do2)
        sys.exit(0)
    vntcall(dirw, ilist, olist, do1, do2, paired, ref_gatk,
            gatk, parallel, temp_dir,
            pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email)


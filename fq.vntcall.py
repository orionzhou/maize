#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import numpy as np
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
    fjob_pre = "50"
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
        pre0 = "%s/%s" % (do1, geno)
        fho1.write("%s -T RealignerTargetCreator -nt %s \
                -R %s %s -o %s.list\n" % \
                (jgatk, '24', ref_gatk, bamstr, pre0))
        fho2.write("%s -T IndelRealigner -R %s \
                %s -targetIntervals %s.list -o %s.bam\n" % \
                (jgatk, ref_gatk, bamstr, pre0, pre0))
        pre1 = "%s/%s.1" % (do2, geno)
        #fho3.write("%s -T UnifiedGenotyper -nt 8 -nct 3 -R %s \
        #        -I %s.bam -o %s.vcf\n" % \
        #        (jgatk, ref_gatk, pre1, pre2))
        fho3.write("%s -T HaplotypeCaller -nct %s -R %s \
                -I %s.bam -o %s.vcf\n" % \
                (jgatk, '24', ref_gatk, pre0, pre1))
        pre2 = "%s/%s.2" % (do2, geno)
        fho4.write("genomeCoverageBed -ibam %s.bam \
                -g $genome/Zmays_v4/15.sizes > %s.cov.tsv\n" % (pre0, pre0))
        fho4.write("%s -T SelectVariants -R %s \
                -V %s.vcf -selectType SNP -o %s.snp.vcf\n" % \
                (jgatk, ref_gatk, pre1, pre2))
        fho4.write("%s -T SelectVariants -R %s \
                -V %s.vcf -selectType INDEL -o %s.indel.vcf\n" % \
                (jgatk, ref_gatk, pre1, pre2))
        pre3 = "%s/%s.3" % (do2, geno)
        fho4.write("%s -T VariantFiltration --logging_level ERROR -R %s \
                -V %s.snp.vcf --filterExpression 'QD < 2.0 \
                || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 \
                || ReadPosRankSum < -8.0 || SOR > 4.0' \
                --filterName 'basic_snp_filter' -o %s.snp.vcf\n" % \
                (jgatk, ref_gatk, pre2, pre3))
        fho4.write("%s -T VariantFiltration --logging_level ERROR -R %s \
                -V %s.indel.vcf --filterExpression 'QD < 2.0 \
                || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
                --filterName 'basic_indel_filter' -o %s.indel.vcf\n" % \
                (jgatk, ref_gatk, pre2, pre3))
        pre4 = "%s/%s.4" % (do2, geno)
        fho4.write("%s -T CombineVariants \
                -R %s --assumeIdenticalSamples \
                -V %s.snp.vcf -V %s.indel.vcf \
                -o %s.vcf\n" \
                % (jgatk, ref_gatk, pre3, pre3, pre4))
        pre5 = "%s/%s.5.het" % (do2, geno)
        fho4.write("#%s -T VariantFiltration --logging_level ERROR -R %s \
                -V %s.vcf --genotypeFilterExpression \
                'isHet == 1 || isHomRef == 1' \
                --genotypeFilterName 'het' -o %s.vcf\n" % \
                (jgatk, ref_gatk, pre4, pre5))
        pre6 = "%s/%s.6.maxdp" % (do2, geno)
        fho4.write("#%s -T VariantFiltration --logging_level ERROR -R %s \
                -V %s.vcf --filterExpression 'DP > 9999' \
                --filterName 'max_dp_filter' -o %s.vcf\n" % \
                (jgatk, ref_gatk, pre5, pre6))
        fho4.write("vcf2tsv.py %s.vcf %s.tsv\n" % (pre4, pre4))

    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    
    pbs_walltimes = pbs_walltime.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_queues = pbs_queue.split(",")
    cmds = [[
    ##cmds.append("module load gatk/3.7.0")
    #cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo4))
        "module load picard/2.3.0",
        "module load java/jdk1.8.0_45",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "bash %s" % fo1
    ], [
        "module load picard/2.3.0",
        "module load java/jdk1.8.0_45",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "bash %s" % fo2
    ], [
        "module load picard/2.3.0",
        "module load java/jdk1.8.0_45",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "bash %s" % fo3
    ], [
        "module load picard/2.3.0",
        "module load java/jdk1.8.0_45",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "bash %s" % fo4
    ]]
    njob = len(cmds)
    assert len(pbs_walltimes) == njob, "not %d jobs" % njob
    assert len(pbs_ppns) == njob, "not %d jobs" % njob

    fjobs = ["%s.%s.pbs" % (fjob_pre, chr(97+i)) for i in range(njob)]
    for i in range(njob):
        temdict = {
                "queue": pbs_queues[i],
                "walltime": pbs_walltimes[i],
                "ppn": pbs_ppns[i],
                "email": pbs_email,
                "cmds": "\n".join(cmds[i])
        }
        fho = open(fjobs[i], "w")
        fho.write(src.substitute(temdict))

    init()
    print(Fore.GREEN)
    print("%s job scripts has been generated: %s" % (njob, ", ".join(fjobs)))
    print("Please check, make necessary changes, then type:")
    print(Fore.RED + "qsub -W depend=afterany:??? %s" % fjobs[1])
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
    import argparse
    import configparser
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Variant calling'
    )
    parser.add_argument(
            'config', nargs = '?', default = "config.ini",
            help = 'config file'
    )
    parser.add_argument(
            '--check', action = "store_true",
            help = 'run the script in check mode'
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


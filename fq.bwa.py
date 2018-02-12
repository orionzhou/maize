#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import configparser
from string import Template
from colorama import init, Fore, Back, Style
import pysam
from crimson import picard

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

def run_bwa(dirw, ilist, olist, do1, do2, do3, paired, db_bwa, ref_gatk,
        bwa, snpEff, samtools, parallel, temp_dir,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2, fo3, fo4, fo5 = "40.1.bwa.sh", "40.2.picard.sh", \
            "40.3.dedup.sh", "40.4.idx.sh", "40.5.stat.sh"
    fho1, fho2, fho3, fho4, fho5 = \
            [open(x, "w") for x in (fo1, fo2, fo3, fo4, fo5)]
    for diro in [do1, do2, do3]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        spec, tiss, geno = row[1:4]
        fsam = "%s/%s.sam" % (do1, sid)
        exist_fsam = check_sam(fsam)
        sam = "_".join(list(filter(None, [spec, geno, tiss])))
        rg = "'@RG\\tID:%s\\tLB:%s\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:%s'" % \
                (sid, sid, sam)
        pre = "%s/%s" % (do1, sid)
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            fho1.write("%s mem -t %s -M -R %s %s %s %s > %s.PE.sam\n" % \
                    (bwa, pbs_ppn, rg, db_bwa, f1p, f2p, pre))
            fho1.write("%s mem -t %s -M -R %s %s %s > %s.1.SE.sam\n" % \
                    (bwa, pbs_ppn, rg, db_bwa, f1u, pre))
            fho1.write("%s mem -t %s -M -R %s %s %s > %s.2.SE.sam\n" % \
                    (bwa, pbs_ppn, rg, db_bwa, f2u, pre))
            fho2.write("$PTOOL/picard.jar MergeSamFiles I=%s.PE.sam \
                    I=%s.1.SE.sam I=%s.2.SE.sam O=%s.bam\n" % \
                    (pre, pre, pre, pre))
        else:
            fr, rc, ft, rrc = row[5:9]
            fho1.write("%s mem -t %s -M -R %s %s %s > %s.sam\n" % \
                    (bwa, pbs_ppn, rg, db_bwa, ft, pre))
            fho2.write("$PTOOL/picard.jar MergeSamFiles \
                    I=%s.sam O=%s.bam\n" % \
                    (pre, pre))
        pre2= "%s/%s" % (do2, sid)
        fho3.write("$PTOOL/picard.jar MarkDuplicates INPUT=%s.bam \
                OUTPUT=%s.bam METRICS_FILE=%s.txt\n" % (pre, pre2, pre2))
        fho4.write("$PTOOL/picard.jar BuildBamIndex INPUT=%s.bam\n" % pre2)
        pre3= "%s/%s" % (do3, sid)
        fho5.write("$PTOOL/picard.jar CollectAlignmentSummaryMetrics \
                R=%s I=%s.bam O=%s.sum.txt\n" % \
                (ref_gatk, pre2, pre3))
        fho5.write("$PTOOL/picard.jar CollectInsertSizeMetrics \
                INPUT=%s.bam OUTPUT=%s.ins.txt HISTOGRAM_FILE=%s.hist.pdf\n" \
                % (pre2, pre3, pre3))
        #fho5.write("%s depth -a %s.bam > %s.depth.txt\n" % (samtools, pre2, pre3))
    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("module load picard/2.3.0")
    cmds.append("export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir)
    cmds.append("bash %s" % fo1)
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo2))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo3))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo4))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo5))
    cmd = "\n".join(cmds)
    
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": cmd
    }
    fo = "40.pbs"
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
def bwa_check(dirw, ilist, olist, do1, do2, do3):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + ["BAM", "Pair_Dup", "Pair_Opt_Dup", \
            "Pct_Dup", \
            "Pair", "Pair_Map", "Pair_Orphan", "Pair_Unmap", \
            "Pair_Map_Hq", "Unpair", "Unpair_Map", "Unpair_Map_Hq"])+"\n")
    else:
        fho.write("\t".join(cols + ["BAM"]) + "\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        bam = "%s/%s.bam" % (do2, sid)
        assert check_bam(bam), "%s not there" % bam
        fd = "%s/%s.txt" % (do2, sid)
        rd = picard.parse(fd)['metrics']['contents']
        fs = "%s/%s.sum.txt" % (do3, sid)
        rs1 = picard.parse(fs)['metrics']['contents']
        rs = { rs1[i]['CATEGORY']: rs1[i] for i in list(range(len(rs1))) }
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            pair_map = rd['READ_PAIRS_EXAMINED']
            pair_dup = rd['READ_PAIR_DUPLICATES']
            pair_opt = rd['READ_PAIR_OPTICAL_DUPLICATES']
            pct_dup = rd['PERCENT_DUPLICATION']
            pair = rs['FIRST_OF_PAIR']['TOTAL_READS']
            assert pair_map == rs['FIRST_OF_PAIR']['READS_ALIGNED_IN_PAIRS'], "error 0"
            pair_map1 = rs['FIRST_OF_PAIR']['PF_READS_ALIGNED']
            pair_map_hq1 = rs['FIRST_OF_PAIR']['PF_HQ_ALIGNED_READS']
            pair_map2 = rs['SECOND_OF_PAIR']['PF_READS_ALIGNED']
            pair_map_hq2 = rs['SECOND_OF_PAIR']['PF_HQ_ALIGNED_READS']
            unpair = rs['UNPAIRED']['TOTAL_READS']
            unpair_map = rs['UNPAIRED']['PF_READS_ALIGNED']
            unpair_map_hq = rs['UNPAIRED']['PF_HQ_ALIGNED_READS']
            pair_orphan = pair_map1 + pair_map2 - pair_map * 2
            pair_unmap = pair - pair_map - pair_orphan
            pair_map_hq = int((pair_map_hq1+pair_map_hq2)/2)
            assert pair == int(rrc), "error 1"
            assert int(rc1)+int(rc2) == unpair, "error 2"
            stats = map(str, [pair_dup, pair_opt, pct_dup, \
                    pair, pair_map, pair_orphan, pair_unmap, \
                    pair_map_hq, unpair, unpair_map, unpair_map_hq])
            fho.write("\t".join(row + [bam] + list(stats)) + "\n")
        else:
            fr, rc, ft, rrc = row[5:9]
            fho.write("\t".join(row + [bam]) + "\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Run BWA-MEM, dedup and collect stats'
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
    cfg = cfg['bwa']
    dirw, ilist, olist, do1, do2, do3 = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['outdir1'], \
            cfg['outdir2'], cfg['outdir3']
    paired = cfg.getboolean('paired')
    annotation = cfg['annotation']
    db_bwa = cfg['db_bwa']
    ref_gatk = cfg['ref_gatk']
    temp_dir = cfg['temp_dir']
    bwa, snpEff, samtools, parallel = \
            cfg['bwa'], cfg['snpEff'], cfg['samtools'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if args.check:
        bwa_check(dirw, ilist, olist, do1, do2, do3)
        sys.exit(0)
    run_bwa(dirw, ilist, olist, do1, do2, do3, paired, db_bwa, ref_gatk,
            bwa, snpEff, samtools, parallel, temp_dir,
            pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email)


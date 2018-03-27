#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import numpy as np
import sys
import logging

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob

from fadapa import Fadapa

def parse_fastqc(fqc):
    assert op.isfile(fqc), "%s not exist" % fqc
    qc = Fadapa(fqc)
    r = dict()
    for ary in qc.clean_data('Basic Statistics'):
        r[ary[0]] = ary[1]
    return r

def fq_trim(cfg, check):
    cfg = cfg['fastq_trim']
    dirw, ilist, olist, jobpre, do1, do2, do3 = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['job_prefix'], \
            cfg['outdir1'], cfg['outdir2'], cfg['outdir3']
    paired = cfg.getboolean('paired')
    temp_dir = cfg['temp_dir']
    f_adp, fastqc, trimm, parallel = \
            cfg['adapter'], cfg['fastqc'], cfg['trimmomatic'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if check:
        fq_trim_check(dirw, ilist, olist, do1, do2, do3, paired)
        sys.exit(0)
    
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(f_adp), "%s not exist" % f_adp 
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2, fo3 = ["%s.%d.sh" % (jobpre, i) for i in range(1,4)]
    fho1, fho2, fho3 = [open(fo, "w") for fo in [fo1, fo2, fo3]]
    assert op.isfile(trimm), "%s is not there" % trimm
    for diro in [do1, do2, do3]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        if paired:
            f1, f2 = row[1:3]
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            f11 = "%s/%s_1.PE.fq.gz" % (do2, sid)
            f12 = "%s/%s_1.SE.fq.gz" % (do2, sid)
            f21 = "%s/%s_2.PE.fq.gz" % (do2, sid)
            f22 = "%s/%s_2.SE.fq.gz" % (do2, sid)
            fho1.write("%s -o %s --extract -f fastq %s %s\n" % \
                    (fastqc, do1, f1, f2))
            fho2.write("java -Xmx2500M -jar %s PE -threads 4 \
                    %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10:8:no \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35\n" % \
                    (trimm, f1, f2, f11, f12, f21, f22, f_adp))
            fho3.write("%s -o %s --extract -f fastq %s %s %s %s\n" % \
                    (fastqc, do3, f11, f12, f21, f22))
        else:
            f1 = row[1]
            assert op.isfile(f1), "%s not there" % f1
            fho1.write("%s -o %s --extract -f fastq %s\n" % \
                    (fastqc, do1, f1))
            fo = "%s/%s.fq.gz" % (do2, sid)
            #fob = "%s/%s_1.fastq.gz" % (do2, sid)
            #if op.isfile(fob):
            #    os.system("mv %s %s" % (fob, fo))
            fho2.write("java -Xmx2500M -jar %s SE -threads 4 \
                    %s %s ILLUMINACLIP:%s:2:30:10:8:no \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35\n" % \
                    (trimm, f1, fo, f_adp))
            fho3.write("%s -o %s --extract -f fastq %s\n" % \
                    (fastqc, do3, fo))
    
    cmds = []
    cmds.append("export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir)
    cmds.append("cd %s" % dirw)
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo1))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo2))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo3))
    
    pbsjob = PbsJob(queue = pbs_queue, 
            ppn = pbs_ppn, 
            walltime = pbs_walltime,
            email = pbs_email,
            cmds = "\n".join(cmds)
    )
    fjob = "%s.pbs" % jobpre
    pbsjob.write(fjob)
    logging.debug("Job script '%s' has been created" % fjob)

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
            f1, f2 = row[1:3]
            f1p = "%s_1.PE.fq.gz" % (sid)
            f1s = "%s_1.SE.fq.gz" % (sid)
            f2p = "%s_2.PE.fq.gz" % (sid)
            f2s = "%s_2.SE.fq.gz" % (sid)
            p1p, p1s, p2p, p2s = ["%s/%s" % (do2, x) for x in [f1p, f1s, f2p, f2s]]
            assert op.isfile(p1p), "%s not exist" % p1p
            assert op.isfile(p1s), "%s not exist" % p1s
            assert op.isfile(p2p), "%s not exist" % p2p
            assert op.isfile(p2s), "%s not exist" % p2s
            pre1 = op.basename(f1).split(".")[0]
            pre2 = op.basename(f2).split(".")[0]
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
            f1 = row[1]
            p1 = "%s/%s.fq.gz" % (do2, sid)
            assert op.isfile(p1), "%s not exist" % p1
            pre1 = op.basename(f1).split(".")[0]
            q1 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre1)
            q2 = "%s/%s_fastqc/fastqc_data.txt" % (do3, sid)
            r1, r2 = parse_fastqc(q1), parse_fastqc(q2)
            rc1, rc2 = r1['Total Sequences'], r2['Total Sequences']
            fho.write("\t".join(row + [rc1, p1, rc2]) + "\n")

def hisat(cfg, check):
    cfg = cfg['hisat']
    dirw, ilist, olist, jobpre, diro1, diro2 = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['job_prefix'], \
            cfg['outdir1'], cfg['outdir2']
    paired = cfg.getboolean('paired')
    temp_dir = cfg['temp_dir']
    ref_gatk = cfg['ref_gatk']
    gatk = cfg['gatk']
    db_hisat, hisat, samtools, parallel = \
            cfg['db_hisat'], cfg['hisat'], cfg['samtools'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if check:
        hisat_check(dirw, ilist, olist, diro1, diro2, paired)
        sys.exit(0)
    
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2, fo2b, fo3 = ["%s.%s.sh" % (jobpre, i) for i in \
            ['1.hisat','2.bam','2.bamidx','3.stat']]
    fho1, fho2, fho2b, fho3 = [open(x, "w") for x in [fo1, fo2, fo2b, fo3]]
    for diro in [diro1, diro2]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    jgatk = "java -jar %s" % gatk
    pbs_queues = pbs_queue.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_walltimes = pbs_walltime.split(",")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        pre1= "%s/%s" % (diro1, sid)
        fsam = "%s.sam" % pre1
        fbam = "%s.bam" % pre1
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[1:11]
            if not op.isfile(fsam):
                fho1.write("%s -p %s -x %s -q -1 %s -2 %s -U %s,%s \
                        --rg-id %s --rg SM:%s -S %s.sam\n" % \
                        (hisat, pbs_ppns[0], db_hisat, f1p, f2p, f1u, f2u, sid, sid, pre1))
        else:
            fr, rc, ft, rrc = row[1:5]
            if not op.isfile(fsam):
                fho1.write("%s -p %s -x %s -q -U %s \
                        --rg-id %s --rg SM:%s -S %s.sam\n" % \
                        (hisat, pbs_ppns[0], db_hisat, ft, sid, sid, pre1))
        if not op.isfile(fbam):
            #fho2.write("$PTOOL/picard.jar SortSam I=%s.sam \
            #        O=%s.bam SORT_ORDER=coordinate\n" % (pre1, pre1))
            #fho2.write("$PTOOL/picard.jar BuildBamIndex INPUT=%s.bam\n" \
            #        % pre1)
            fho2.write("%s sort -m 2500M -O bam -o %s.bam %s.sam\n" % (samtools, pre1, pre1))
            fho2b.write("%s index %s.bam\n" % (samtools, pre1))
        pre2 = "%s/%s" % (diro2, sid)
        #fho3.write("%s -T IndelRealigner -R %s \
        #        -I %s.bam -U ALLOW_N_CIGAR_READS \
        #        -targetIntervals %s -known %s -o %s.bam\n" % \
        #        (jgatk, ref_gatk, pre1, frta, fvcf, pre2))
        fho3.write("$PTOOL/picard.jar CollectAlignmentSummaryMetrics \
                R=%s I=%s.bam O=%s.sum.txt\n" % \
                (ref_gatk, pre1, pre2))
        fho3.write("$PTOOL/picard.jar CollectInsertSizeMetrics \
                INPUT=%s.bam OUTPUT=%s.ins.txt HISTOGRAM_FILE=%s.hist.pdf\n" \
                % (pre1, pre2, pre2))
    
    cmds = [[
        "cd %s" % dirw,
        "bash %s" % fo1
    ], [
        "cd %s" % dirw,
        "%s -j %s < %s" % (parallel, pbs_ppns[1], fo2),
        "%s -j %s < %s" % (parallel, pbs_ppns[1], fo2b)
    ], [
        "module load picard/2.3.0",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "%s -j %s < %s" % (parallel, pbs_ppns[2], fo3)
    ]]
    njob = len(cmds)
    assert len(pbs_walltimes) == njob, "not %d jobs" % njob
    assert len(pbs_ppns) == njob, "not %d jobs" % njob

    fjobs = ["%s.%s.pbs" % (jobpre, chr(97+i)) for i in range(njob)]
    for i in range(njob):
        pbsjob = PbsJob(queue = pbs_queues[i],
                ppn = pbs_ppns[i],
                walltime = pbs_walltimes[i],
                email = pbs_email,
                cmds = "\n".join(cmds[i])
        )
        fjob = "%s.pbs" % jobpre
        pbsjob.write(fjobs[i])
        
    logging.debug("%s job scripts were created: %s" % (njob, ", ".join(fjobs)))
    logging.debug("qsub %s" % fjobs[0])
    logging.debug("qsub -W depend=afterok:??? %s" % fjobs[1])

def hisat_check(dirw, ilist, olist, diro1, diro2, paired):
    from crimson import picard
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + ["BAM",
            "Pair", "Pair_Map", "Pair_Orphan", "Pair_Unmap", \
            "Pair_Map_Hq", "Unpair", "Unpair_Map", "Unpair_Map_Hq"])+"\n")
    else:
        fho.write("\t".join(cols + ["BAM",
            "Total", "Mapped", "Mapped_Hq"]) + "\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        bam = "%s/%s.bam" % (diro1, sid)
        assert op.isfile(bam), "%s not exist" % bam
        fs = "%s/%s.sum.txt" % (diro2, sid)
        rs1 = picard.parse(fs)['metrics']['contents']
        if type(rs1) == dict: rs1 = [rs1]
        rs = { rs1[i]['CATEGORY']: rs1[i] for i in list(range(len(rs1))) }
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[1:11]
            pair = rs['FIRST_OF_PAIR']['TOTAL_READS']
            pair_map = rs['FIRST_OF_PAIR']['READS_ALIGNED_IN_PAIRS']
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
            assert pair == int(rrc), "error 1: %d %s" % (pair, rrc)
            assert int(rc1)+int(rc2) == unpair, "error 2"
            stats = map(str, [pair, pair_map, pair_orphan, pair_unmap, \
                    pair_map_hq, unpair, unpair_map, unpair_map_hq])
            fho.write("\t".join(row + [bam] + list(stats)) + "\n")
        else:
            fr, rc, ft, rrc = row[1:5]
            unpair = rs['UNPAIRED']['TOTAL_READS']
            unpair_map = rs['UNPAIRED']['PF_READS_ALIGNED']
            unpair_map_hq = rs['UNPAIRED']['PF_HQ_ALIGNED_READS']
            stats = map(str, [unpair, unpair_map, unpair_map_hq])
            fho.write("\t".join(row + [bam] + list(stats)) + "\n")

def htseq(cfg, check):
    cfg = cfg['htseq']
    dirw, ilist, olist, jobpre, diro = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['job_prefix'], \
            cfg['outdir']
    paired = cfg.getboolean('paired')
    srd, annotation = cfg['stranded'], cfg['annotation']
    htseq, samtools, parallel = cfg['htseq'], cfg['samtools'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if check:
        htseq_check(dirw, ilist, olist, diro)
        sys.exit(0)
    
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1 = "%s.1.htseq.sh" % jobpre
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
            fbam = row[11]
            if not op.isfile(fsen) or os.stat(fsen).st_size == 0:
                fho1.write("%s view -f 1 -F 256 %s | %s -r pos -s %s \
                        -t exon -i gene_id -m union -a 20 - %s > %s/%s.txt\n" % 
                        (samtools, fbam, htseq, 'reverse', annotation, diro, sid))
            if not op.isfile(fant) or os.stat(fant).st_size == 0:
                fho1.write("%s view -f 1 -F 256 %s | %s -r pos -s %s \
                        -t exon -i gene_id -m union -a 20 - %s > %s/%s.as.txt\n" % 
                        (samtools, fbam, htseq, 'yes', annotation, diro, sid))
        else:
            fbam = row[5]
            if not op.isfile(fsen) or os.stat(fsen).st_size == 0:
                fho1.write("%s view -F 256 %s | %s -r pos -s no \
                        -t exon -i gene_id -m union -a 20 - %s > %s/%s.txt\n" % 
                        (samtools, fbam, htseq, annotation, diro, sid))

    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("module load python2")
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo1))
    cmd = "\n".join(cmds)

    pbsjob = PbsJob(queue = pbs_queue, 
            ppn = pbs_ppn, 
            walltime = pbs_walltime,
            email = pbs_email,
            cmds = "\n".join(cmds)
    )
    fjob = "%s.pbs" % jobpre
    pbsjob.write(fjob)
    logging.debug("Job script '%s' has been created" % fjob)

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

def run_ase(cfg, check):
    import pysam
    cfg = cfg['ase']
    dirw, ilist, olist, jobpre, diro = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['job_prefix'], \
            cfg['outdir']
    f_fas = cfg['genome']
    paired = cfg.getboolean('paired')
    samtools, bcftools, parallel = \
            cfg['samtools'], cfg['bcftools'], cfg['parallel']
    target_vcf, gene_bed = cfg['targetvcf'], cfg['gene_bed']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if check:
        ase_check(dirw, ilist, olist, diro, paired)
        sys.exit(0)

    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    dirj = "%s.jobs" % jobpre
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
            fbam = row[11]
        else:
            fbam = row[5]
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
            "bed.ase.sum.py %s.5.tsv %s.8.bed %s.tsv" % (pre, pre, pre),
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

    cmds = []
    cmds.append("cd %s" % dirw),
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fj))
    
    pbsjob = PbsJob(queue = pbs_queue, 
            ppn = pbs_ppn, 
            walltime = pbs_walltime,
            email = pbs_email,
            cmds = "\n".join(cmds)
    )
    fjob = "%s.pbs" % jobpre
    pbsjob.write(fjob)
    logging.debug("Job script '%s' has been created" % fjob)

def run_ase1(dirw, ilist, olist, diro, paired, f_fas, ref_gatk,
        gatk, samtools, parallel, temp_dir,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2 = "41.1.split.sh", "41.2.hc.sh"
    fjob_pre = "41"
    fho1, fho2 = [open(x, "w") for x in [fo1, fo2]]
    for diro in [diro]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    bams = []
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        if paired:
            fbam = row[15]
        else:
            fbam = row[9]
    #    bams.append("-I %s" % fbam)
        bams.append("I=%s" % fbam)
    jgatk = "java -jar %s" % gatk #-Xmx62g
    bamstr = " ".join(bams)
    pre1 = "%s/11" % diro
    fho1.write("$PTOOL/picard.jar MergeSamFiles %s \
            O=%s.bam CREATE_INDEX=true\n" % \
            (bamstr, pre1))
    ##fho1.write("$PTOOL/picard.jar MarkDuplicates %s \
    ##        O=%s.bam M=%s.txt CREATE_INDEX=true\n" % \
    ##        (bamstr, pre1, pre1))
    pre4 = "%s/14.splitn" % diro
    fho1.write("%s -T SplitNCigarReads \
            -R %s -I %s.bam -o %s.bam -rf ReassignOneMappingQuality \
            -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS\n" % \
            (jgatk, ref_gatk, pre1, pre4))
    fho2.write("%s -T HaplotypeCaller -nct 24 \
            -R %s -dontUseSoftClippedBases \
            --maxReadsInRegionPerSample 10000 \
            -stand_call_conf 20 -I %s.bam -o %s/21.raw.vcf\n" % \
            (jgatk, ref_gatk, pre4, diro))
    #fho1.write("%s mpileup -ugf %s %s | bcftools call -vmO v -o %s" \
    #        % (samtools, f_fas, bamstr, fv))
    #fho1.write("%s -T RealignerTargetCreator -nt %s \
    #        --filter_reads_with_N_cigar \
    #        -R %s %s -o %s/01.rt.list\n" % \
    #        (jgatk, pbs_ppn, ref_gatk, bamstr, diro))
    #fho1.write("%s -T IndelRealigner -R %s \
    #        --filter_reads_with_N_cigar \
    #        %s -targetIntervals %s/01.rt.list -o %s/02.realigned.bam\n" % \
    #        (jgatk, ref_gatk, bamstr, diro, diro))
    #fho1.write("%s -T UnifiedGenotyper -nt 8 -nct 3 -R %s \
    #        --filter_reads_with_N_cigar \
    #        %s/02.realigned.bam -o %s/11.raw.vcf\n" % \
    #        (jgatk, ref_gatk, diro, diro))
    #fho1.write("vcf2tsv.py 11.raw.vcf 11.raw.tsv\n")
    
    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    
    pbs_walltimes = pbs_walltime.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_queues = pbs_queue.split(",")
    cmds = [[
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
def run_ase2(dirw, ilist, olist, diro, paired, f_fas, target_vcf,
        samtools, bcftools, parallel,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1 = "24.1.bcftools.sh"
    fjob_pre = "24"
    fho1 = open(fo1, "w")
    #fho1 = [open(x, "w") for x in [fo1]]
    for diro in [diro]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        #if sid > 'BR041':
        #    continue
        if paired:
            fbam = row[15]
        else:
            fbam = row[9]
        fho1.write("%s mpileup -Ou -f %s %s | \
                %s call -C alleles -m -T %s -O v | vcf2ase.py - %s/%s.tsv\n" % 
                (bcftools, f_fas, fbam, bcftools, target_vcf, diro, sid))
    
    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    
    pbs_walltimes = pbs_walltime.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_queues = pbs_queue.split(",")
    cmds = [[
        "cd %s" % dirw,
        "%s -j %s < %s" % (parallel, pbs_ppns[0], fo1)
    ]
    ]
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
    #print(Fore.RED + "qsub -W depend=afterany:??? %s" % fjobs[1])
    print(Style.RESET_ALL)
def ase_check(dirw, ilist, olist, diro, paired):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + ["BAM",
            "Pair", "Pair_Map", "Pair_Orphan", "Pair_Unmap", \
            "Pair_Map_Hq", "Unpair", "Unpair_Map", "Unpair_Map_Hq"])+"\n")
    else:
        fho.write("\t".join(cols + ["BAM"]) + "\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        bam = "%s/%s.bam" % (diro, sid)
        assert check_bam(bam), "%s not exist" % bam
        fs = "%s/%s.sum.txt" % (diro, sid)
        rs1 = picard.parse(fs)['metrics']['contents']
        rs = { rs1[i]['CATEGORY']: rs1[i] for i in list(range(len(rs1))) }
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            pair = rs['FIRST_OF_PAIR']['TOTAL_READS']
            pair_map = rs['FIRST_OF_PAIR']['READS_ALIGNED_IN_PAIRS']
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
            stats = map(str, [pair, pair_map, pair_orphan, pair_unmap, \
                    pair_map_hq, unpair, unpair_map, unpair_map_hq])
            fho.write("\t".join(row + [bam] + list(stats)) + "\n")
        else:
            fr, rc, ft, rrc = row[5:9]
            fho.write("\t".join(row + [bam]) + "\n")

if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Illumina RNA-Seq pipeline'
    )
    parser.add_argument('config', nargs = "?", default = "config.ini", help = 'config file')
    parser.add_argument('--check', action = "store_true", help = 'run script in check mode')
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("fqtrim", help = "Trimming and QC fastq files")
    sp1.set_defaults(func = fq_trim)
    
    sp1 = sp.add_parser("hisat", help = 'Map fastq seqs to genome using hisat2')
    sp1.set_defaults(func = hisat)

    sp1 = sp.add_parser("htseq", help = 'Quantify gene expression using htseq')
    sp1.set_defaults(func = htseq)
    
    sp1 = sp.add_parser("ase", help = 'Allele-specific Expression calculation')
    sp1.set_defaults(func = run_ase)
    
    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    if args.command:
        args.func(cfg, args.check)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


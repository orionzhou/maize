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


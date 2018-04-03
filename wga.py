#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import pysam

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open, ndigit
from maize.formats.pbs import PbsJob

def prepare(cfg):
    cfg = cfg['prepare']
    dirw, qry, tgt = cfg['dirw'], cfg['qry'], cfg['tgt']
    qry_fas, tgt_fas = cfg['qry_fas'], cfg['tgt_fas']
    tmpdir = cfg['temp_dir']
    npieces = int(cfg['npieces'])
    diro1 = cfg['outdir1']
    
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

    sh("genome %s/%s fasta" % (dirw, "01_tgt_genome"))
    sh("genome %s/%s blat" % (dirw, "01_tgt_genome"))
    sh("genome %s/%s fasta" % (dirw, "02_qry_genome"))
    sh("genome %s/%s blat" % (dirw, "02_qry_genome"))
    
    sh("bed filter --minsize 1000 02_qry_genome/16.gap.bed > 04.qry.gap.bed")
    sh("subtractBed -nonamecheck -a 02_qry_genome/15.bed -b 04.qry.gap.bed | bed filter -min 100 - | bed makewindow -w 100000 -s 95000 - > 05.qry.clean.bed")
    sh("bed size 05.qry.clean.bed")
    sh("fasta extract 02_qry_genome/11_genome.fas 05.qry.clean.bed > 06.qry.fas")
    sh("fasta split --N %d %s %s" % (npieces, '06.qry.fas', diro1))

def run_blat(cfg):
    cfg = cfg['blat']
    dirw, jobpre, diro1, diro2 = \
            cfg['dirw'], cfg['job_prefix'], cfg['outdir1'], cfg['outdir2']
    qry, tgt = cfg['qry'], cfg['tgt']
    parallel = cfg['parallel']
    npieces, npieces2 = int(cfg['npieces']), int(cfg['npieces2'])
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
   
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
    for subdir in [diro1, diro2]:
        if not op.isdir(subdir):
            mkdir(subdir)
    
    dirt = "01_tgt_genome"
    dirq = "02_qry_genome"
    tgt_fas = "%s/11_genome.fas" % dirt
    qry_fas = "%s/11_genome.fas" % dirq
    tgt_2bit = "%s/21.blat/db.2bit" % dirt
    qry_2bit = "%s/21.blat/db.2bit" % dirq
    tgt_ooc = "%s/21.blat/db.2bit.tile11.ooc" % dirt
    tgt_size = "%s/15.sizes" % dirt
    qry_size = "%s/15.sizes" % dirq
    tgt_size_bed = "%s/15.bed" % dirt
    qry_size_bed = "%s/15.bed" % dirq
    tgt_gap = "%s/16.gap.bed" % dirt
    qry_gap = "%s/16.gap.bed" % dirq
    for fn in [tgt_fas, qry_fas, tgt_2bit, qry_2bit, tgt_ooc,
            tgt_size, qry_size, tgt_size_bed, qry_size_bed, tgt_gap, qry_gap]:
        if not op.isfile(fn):
            logging.error("%s is not there" % fn)
            sys.exit()
    
    pbs_queues = pbs_queue.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_walltimes = pbs_walltime.split(",")
    njob = len(pbs_queues)
    assert len(pbs_walltimes) == njob, "not %d jobs" % njob
    assert len(pbs_ppns) == njob, "not %d jobs" % njob
    fbs = ["%sb.%d.sh" % (jobpre, i+1) for i in range(njob)]
    fjs = ["%sj.%d.pbs" % (jobpre, i+1) for i in range(njob)]
    bcmds, jcmds = [], []

    #1 blat
    cmds = []
    bcmds.append(cmds)
    prepre = "part.%%0%dd" % ndigit(npieces-1)
    jcmds.append([
        "let i=${PBS_ARRAYID}",
        "cd %s" % dirw,
        "printf -v pre %s \"$i\"" % prepre,
        "pblat %s %s/${pre}.fas -threads=%s -ooc=%s %s/${pre}.psl" % \
            (tgt_2bit, diro1, pbs_ppns[0], tgt_ooc, diro2)
    ])

    #2 process blat 
    bcmds.append([
        "pslCat -nohead 12.blat/part.*.psl > 12.psl",
        "psl qcoord 12.psl %s > 13.psl" % qry_size,
        "pslCheck -querySizes=%s -targetSizes=%s -pass=14.check.psl 13.psl" %
            (qry_size, tgt_size),
        "axtChain -linearGap=medium -psl 14.check.psl %s %s 21.chain" %
            (tgt_2bit, qry_2bit),
        "chainPreNet 21.chain %s %s 23.chain" % (tgt_size, qry_size),
        "chain 2bed --qry 23.chain | sortBed -i stdin | \
                mergeBed -i stdin > 23.bed",
        "subtractBed -a %s -b %s -nonamecheck | \
                subtractBed -a stdin -b 23.bed -nonamecheck | \
                bed filter --minsize 50 - > 24.nov.bed" %
                (qry_size_bed, qry_gap),
        "seqret.pl -d %s -b 24.nov.bed -o 24.nov.fas" % qry_fas,
        "rm 23.chain 23.bed",
        "fasta split --N %d %s %s" % (npieces2, '24.nov.fas', '25.nov'),
    ])
    jcmds.append([
        "cd %s" % dirw,
        "bash %s" % fbs[1],
    ])
   
    #3 blat nov
    cmds = []
    bcmds.append(cmds)
    prepre = "part.%%0%dd" % ndigit(npieces2-1)
    jcmds.append([
        "let i=${PBS_ARRAYID}",
        "cd %s" % dirw,
        "printf -v pre %s \"$i\"" % prepre,
        "pblat %s %s/${pre}.fas -threads=%s -ooc=%s %s/${pre}.psl" % \
            (tgt_2bit, '25.nov', pbs_ppns[2], tgt_ooc, '25.nov')
    ])

    #4 process blat
    bcmds.append([
        "pslCat -nohead 25.nov/part.*.psl > 25.nov.psl",
        "psl qcoord 25.nov.psl %s > 26.psl" % qry_size,
        "pslCheck -querySizes=%s -targetSizes=%s -pass=27.check.psl 26.psl" %
            (qry_size, tgt_size),
        "pslCat 14.check.psl 27.check.psl > 31.1.psl",
        "pslSwap 31.1.psl 41.1.psl",
        "rm 25.nov.psl 26.psl",

        "axtChain -linearGap=medium -psl 31.1.psl %s %s 31.2.chain" % (tgt_2bit, qry_2bit),
        "chainPreNet 31.2.chain %s %s 31.3.chain" % (tgt_size, qry_size),
        "chainSwap 31.3.chain 31.3.q.chain",
        "chainNet 31.3.chain %s %s 31.5.net 31.5.q.net" % (tgt_size, qry_size),
        "netChainSubset 31.5.net 31.3.chain stdout | chainSort stdin 31.5.chain",
        "netChainSubset 31.5.q.net 31.3.q.chain stdout | chainSort stdin 31.5.q.chain",
        "chainNet 31.5.q.chain %s %s /dev/null 31.8.net" % (qry_size, tgt_size),
        "netChainSubset 31.8.net 31.3.chain 31.8.chain",
        
        "axtChain -linearGap=medium -psl 41.1.psl %s %s 41.2.chain" % (qry_2bit, tgt_2bit),
        "chainPreNet 41.2.chain %s %s 41.3.chain" % (qry_size, tgt_size),
        "chainSwap 41.3.chain 41.3.q.chain",
        "chainNet 41.3.chain %s %s 41.5.net 41.5.q.net" % (qry_size, tgt_size),
        "netChainSubset 41.5.net 41.3.chain stdout | chainSort stdin 41.5.chain",
        "netChainSubset 41.5.q.net 41.3.q.chain stdout | chainSort stdin 41.5.q.chain",
        "chainNet 41.5.q.chain %s %s /dev/null 41.8.net" % (tgt_size, qry_size),
        "netChainSubset 41.8.net 41.3.chain 41.8.chain",
    ])
    jcmds.append([
        "cd %s" % dirw,
        "bash %s" % fbs[3],
    ])

    #5 process vnt
    bcmds.append([
        "snp2vcf.pl -i snp -o snp.vcf -s %s" % qry,
    ])
    jcmds.append([
        "cd %s/31.9" % dirw,
        "bash %s" % fbs[4],
    ])
    
    assert len(bcmds) == njob, "not %d jobs" % njob
    assert len(jcmds) == njob, "not %d jobs" % njob
    for i in range(njob):
        fb, fj = fbs[i], fjs[i]
        if op.isfile(fb):
            os.remove(fb)
        if len(bcmds[i]) > 0:
            fhb = open(fb, "w")
            fhb.write("\n".join(bcmds[i]) + "\n")
            fhb.close()
        pbsjob = PbsJob(queue = pbs_queues[i],
                ppn = pbs_ppns[i],
                walltime = pbs_walltimes[i],
                email = pbs_email,
                cmds = "\n".join(jcmds[i])
        )
        pbsjob.write(fj)
        
    logging.debug("%s job scripts were created" % njob)
    eprint("qsub -t 0-%d %s" % (npieces-1, fjs[0]))
    eprint("qsub -W depend=afterok: %s" % fjs[1])
    eprint("qsub -t 0-%d %s" % (npieces2-1, fjs[2]))
    eprint("qsub -W depend=afterok: %s" % fjs[3])
    eprint("qsub -W depend=afterok: %s" % fjs[4])

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
    
    sp1 = sp.add_parser("blat", help = "generate blat jobs")
    sp1.set_defaults(func = run_blat)
    
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


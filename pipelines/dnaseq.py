#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re
import logging
from astropy.table import Table, Column

from maize.apps.base import AttrDict, str2bool, eprint, sh, mkdir, which
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob, create_job_chain

def check_cfg_mapping(c):
    c.outdirs = c.outdir.split(",")
    assert len(c.outdirs) == 2, "not 2 outdirs: %s" % c.outdir
    
    for subdir in [c.dirw, c.temp_dir] + c.outdirs:
        if not op.isdir(subdir):
            mkdir(subdir)
    
    for fn in [c.ilist, c.genome, c.gff]:
        assert op.isfile(fn), "cannot read %s" % fn

    for key in 'samtools parallel sambamba bcftools bedtools'.split():
        fp = which(c[key])
        assert fp is not None, "not executable: %s" % c[key]
        c[key] = fp

    c.paired = str2bool(c.paired)

    if c.mapper == 'bwa':
        c.bwa = which(c.bwa)
        assert c.bwa is not None, "not executable: %s" % c.bwa
    elif c.mapper == 'hisat2':
        c.hisat2 = which(c.hisat2)
        assert c.hisat2 is not None, "not executable: %s" % c.hisat2
    elif c.mapper == 'bowtie2':
        c.bowtie2 = which(c.bowtie2)
        assert c.bowtie2 is not None, "not executable: %s" % c.bowtie2
    else:
        logging.error("unsupported mapper: %s" % c.mapper)
        sys.exit(1)
    
    njob = 3
    c.pbs_walltimes = c.pbs_walltime.split(",")
    c.pbs_ppns = c.pbs_ppn.split(",")
    c.pbs_queues = c.pbs_queue.split(",")
    assert njob == len(c.pbs_queues) == len(c.pbs_walltimes) == len(c.pbs_ppns), "not %d jobs: %s" % (njob, c.pbs_queue)
    c.njob = njob

    return c

def mapping(cfg, args):
    c = AttrDict(cfg['mapping'])
    c = check_cfg_mapping(c)
    if args.check:
        mapping_check(c)
        return 0
    os.chdir(c.dirw)

    jcmds = [[
        "cd %s" % c.dirw,
    ], [
        "cd %s" % c.dirw,
    ], [
        "cd %s" % c.dirw,
    ]]
    bcfgs = [
        [dict(opt = 'bash')], 
        [dict(opt = 'parallel', thread = c.pbs_ppns[1])],
        [dict(opt = 'bash'),
        dict(opt = 'parallel', thread = c.pbs_ppns[2]),
        ],
    ]

    assert c.njob == len(bcfgs) == len(jcmds), "not %d jobs" % c.njob
    jobs = []
    for i in range(c.njob):
        prefix = "%s.%d" % (c.job_prefix, i+1)
        jcfg = {
            'queue': c.pbs_queues[i],
            'ppn': c.pbs_ppns[i], 
            'walltime': c.pbs_walltimes[i],
            'email': c.pbs_email,
        }
        job = PbsJob.from_cfg(jcfg = jcfg, jcmds = jcmds[i], bcfgs = bcfgs[i],
                prefix = prefix, njob = len(bcfgs[i]), 
                bash = c.bash, parallel = c.parallel)
        jobs.append(job)
 
    t = Table.read(c.ilist, format = 'ascii.tab')
    nrow = len(t)
    for i in range(nrow):
        sid = t['sid'][i]
        pre1= "%s/%s" % (c.outdirs[0], sid)
        fsam = "%s.sam" % pre1
        input_str = ''
        if c.paired:
            f1p = t["TrimmedReadFile1Paired"][i]
            f1u = t["TrimmedReadFile1Unpaired"][i]
            f2p = t["TrimmedReadFile2Paired"][i]
            f2u = t["TrimmedReadFile2Unpaired"][i]
            if c.mapper == 'hisat2' or c.mapper == 'bowtie2':
                input_str = "-1 %s -2 %s -U %s,%s" % (f1p, f2p, f1u, f2u)
            elif c.mapper == 'bwa':
                input_str = "%s %s" % (f1p, f2p)
        else:
            ft = t["TrimmedReadFile"][i]
            if c.mapper == 'hisat2' or c.mapper == 'bowtie2':
                input_str = "-U %s" % ft
            elif c.mapper == 'bwa':
                input_str = "%s" % ft
        if c.mapper == 'bwa':
            jobs[0].subjobs[0].add_cmd("%s mem -t %s %s %s \
                    -R '@RG\\tID:%s\\tSM:%s' -a > %s.sam" % \
                    (c.bwa, c.pbs_ppns[0], c.bwa_db, input_str, \
                    sid, sid, pre1))
        elif c.mapper == 'hisat2':
            jobs[0].subjobs[0].add_cmd("%s -p %s -x %s -q %s \
                    --no-spliced-alignment --rg-id %s --rg SM:%s -S %s.sam" % \
                    (c.hisat2, c.pbs_ppns[0], c.hisat_db, input_str, \
                    sid, sid, pre1))
        elif c.mapper == 'bowtie2':
            jobs[0].subjobs[0].add_cmd("%s -p %s -x %s -q %s \
                    --rg-id %s --rg SM:%s --sensitive -S %s.sam" % \
                    (c.bowtie2, c.pbs_ppns[0], c.bowtie_db, input_str, \
                    sid, sid, pre1))
        
        fbam = "%s.bam" % pre1
        jobs[1].subjobs[0].add_cmd("%s view -Sb %s.sam -o %s.raw.bam" % \
                (c.samtools, pre1, pre1))
        jobs[2].subjobs[0].add_cmd("%s sort -t %s -m 60GB %s.raw.bam -o %s.bam" % \
                (c.sambamba, c.pbs_ppns[2], pre1, pre1))
        #bcmds[2].append("%s index -t %s %s.bam" % (sambamba, pbs_ppns[2], pre1))
    
        pre2 = "%s/%s" % (c.outdirs[1], sid)
        jobs[2].subjobs[1].add_cmd("bam stat %s.bam --isize %s.ins.tsv > %s.tsv" % \
                (pre1, pre2, pre2))
  
    for job in jobs:
        job.write()
    fj = "%s.sh" % c.job_prefix
    create_job_chain([job.fname for job in jobs], fj)
    logging.debug("job chain with %s jobs was created: %s" % (c.njob, fj))

def mapping_check(cfg):
    t = Table.read(ilist, format = 'ascii.tab')
    nrow = len(t)
    newcols = ''
    if c.paired:
        newcols = '''BAM Pair Pair_Map Pair_Orphan Pair_Unmap
            Pair_Map_Hq Unpair Unpair_Map Unpair_Map_Hq'''.split()
    else: 
        newcols = '''BAM Total Mapped Mapped_Hq'''.split()
    for newcol in newcols:
        t.add_column(Column(name = newcol, length = nrow, dtype = object))
    
    for i in range(nrow):
        sid = t['sid'][i]
        bam = "%s/%s.bam" % (c.outdirs[0], sid)
        assert op.isfile(bam), "%s not exist" % bam
        fs = "%s/%s.tsv" % (c.outdirs[1], sid)
        assert op.isfile(fs), "%s not exist" % fs
        if c.paired:
            t['BAM'][i] = 0#bam
            t['Pair'][i] = 0#pair
            t['Pair_Map'][i] = 0#pair_map
            t['Pair_Orphan'][i] = 0#pair_orphan
            t['Pair_Unmap'][i] = 0#pair_unmap
            t['Pair_Map_Hq'][i] = 0#pair_map_hq
            t['Unpair'][i] = 0#unpair
            t['Unpair_Map'][i] = 0#unpair_map
            t['Unpair_Map_Hq'][i] = 0#unpair_map_hq
        else:
            t['BAM'] = 0#bam
            t['Total'] = 0#unpair
            t['Mapped'] = 0#unpair_map
            t['Mapped_Hq'] = 0#unpair_map_hq
    t.write(t.olist, format='ascii.tab', overwrite=True)


if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Illumina DNA-Seq pipeline(s)'
    )
    parser.add_argument('--config', "--cfg", default = "config.ini", help = 'config file')
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("mapping",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'mapping'
    )
    sp1.add_argument("--check", action = 'store_true', help = "run script in check mode")
    sp1.set_defaults(func = mapping)

    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    if args.command:
        args.func(cfg, args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


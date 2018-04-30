#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from astropy.table import Table, Column

from maize.apps.base import AttrDict, str2bool, eprint, sh, mkdir, which
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob, create_job_chain

from fadapa import Fadapa

def parse_fastqc(fqc):
    assert op.isfile(fqc), "%s not exist" % fqc
    qc = Fadapa(fqc)
    r = dict()
    for ary in qc.clean_data('Basic Statistics'):
        r[ary[0]] = ary[1]
    return r

def check_cfg_fqtrim(c, njob = 1, noutdir = 3):
    c.outdirs = c.outdir.split(",")
    assert len(c.outdirs) == noutdir, "not %s outdirs: %s" % (noutdir, c.outdir)
    
    for subdir in [c.dirw, c.temp_dir] + c.outdirs:
        if not op.isdir(subdir):
            mkdir(subdir)
    
    for fn in [c.ilist, c.adapter, c.trimmomatic]:
        assert op.isfile(fn), "cannot read %s" % fn

    for key in ['fastqc', 'parallel']:
        fp = which(c[key])
        assert fp is not None, "not executable: %s" % c[key]
        c[key] = fp

    c.paired = str2bool(c.paired)
    
    c.pbs_walltimes = c.pbs_walltime.split(",")
    c.pbs_ppns = c.pbs_ppn.split(",")
    c.pbs_queues = c.pbs_queue.split(",")
    assert njob == len(c.pbs_queues) == len(c.pbs_walltimes) == len(c.pbs_ppns), "not %d jobs: %s" % (njob, c.pbs_queue)
    c.njob = njob

    return c
 
def fq_trim(cfg, args):
    c = AttrDict(cfg['fastq_trim'])
    c = check_cfg_fqtrim(c)
    if args.check:
        fq_trim_check(c)
        return 0
    os.chdir(c.dirw)
 
    jcmds = [[
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % c.temp_dir,
        "cd %s" % c.dirw
    ]]
    bcfgs = [
        [dict(opt = 'parallel', thread = c.pbs_ppns[0]),
        dict(opt = 'parallel', thread = c.pbs_ppns[0]),
        dict(opt = 'parallel', thread = c.pbs_ppns[0])]
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
    if c.paired:
        for i in range(nrow):
            sid, f1, f2 = t['sid'][i], t['Readfile1'][i], t['Readfile2'][i]
            assert op.isfile(f1), "%s not there" % f1
            assert op.isfile(f2), "%s not there" % f2
            jobs[0].subjobs[0].add_cmd("%s -o %s --extract -f fastq %s %s" % \
                    (c.fastqc, c.outdirs[0], f1, f2))
            f11, f12, f21, f22 = ["%s/%s_%s.fq.gz" % (c.outdirs[1], sid, x) \
                    for x in ['1.PE', '1.SE', '2.PE', '2.SE']]
            jobs[0].subjobs[1].add_cmd("java -Xmx2500M -jar %s PE -threads 4 \
                    %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10:8:no \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" % \
                    (c.trimmomatic, f1, f2, f11, f12, f21, f22, c.adapter))
            jobs[0].subjobs[2].add_cmd("%s -o %s --extract -f fastq %s %s %s %s" % \
                    (c.fastqc, c.outdirs[2], f11, f12, f21, f22))
    else:
        for i in range(nrow):
            sid, f1 = t['sid'][i], t['Readfile'][i]
            assert op.isfile(f1), "%s not there" % f1
            jobs[0].subjobs[0].add_cmd("%s -o %s --extract -f fastq %s" % \
                    (c.fastqc, c.outdirs[0], f1))
            fo = "%s/%s.fq.gz" % (c.outdirs[1], sid)
            jobs[0].subjobs[1].add_cmd("java -Xmx2500M -jar %s SE -threads 4 \
                    %s %s ILLUMINACLIP:%s:2:30:10:8:no \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35" % \
                    (c.trimmomatic, f1, fo, c.adapter))
            jobs[0].subjobs[2].add_cmd("%s -o %s --extract -f fastq %s" % \
                    (c.fastqc, c.outdirs[2], fo))
    
    for job in jobs:
        job.write()
    fj = "%s.sh" % c.job_prefix
    create_job_chain([job.fname for job in jobs], fj)
    logging.debug("job chain with %s jobs was created: %s" % (c.njob, fj))

def fq_trim_check(c):
    t = Table.read(c.ilist, format = 'ascii.tab')
    nrow = len(t)
    newcols = ''
    if c.paired:
        newcols = '''ReadPairCount
            TrimmedReadFile1Paired
            TrimmedReadFile1Unpaired
            TrimmedReadFile2Paired
            TrimmedReadFile2Unpaired
            RetainedReadPairCount
            RetainedRead1Count
            RetainedRead2Count'''.split()
    else:
        newcols = "ReadCount TrimmedReadFile TrimmedReadCount".split()
    for newcol in newcols:
        t.add_column(Column(name = newcol, length = nrow, dtype = object))

    if c.paired:
        for i in range(nrow):
            sid, f1, f2 = t['sid'][i], t['Readfile1'][i], t['Readfile2'][i]
            f1p, f1s, f2p, f2s = ["%s/%s_%s.fq.gz" % (c.outdirs[1], sid, x) \
                    for x in ['1.PE', '1.SE', '2.PE', '2.SE']]
            assert op.isfile(f1p) and op.isfile(f1s) and op.isfile(f2p) and op.isfile(f2s), "%s not exist" % sid
            pre1 = op.basename(f1).split(".")[0]
            pre2 = op.basename(f2).split(".")[0]
            q11 = "%s/%s_fastqc/fastqc_data.txt" % (c.outdirs[0], pre1)
            q12 = "%s/%s_fastqc/fastqc_data.txt" % (c.outdirs[0], pre2)
            r11, r12 = parse_fastqc(q11), parse_fastqc(q12)
            rc11, rc12 = r11['Total Sequences'], r12['Total Sequences']
            assert rc11 == rc12, "%s: read1 %d != read2 %d" % \
                    (sid, rc11, rc12)
            q21 = "%s/%s_1.PE_fastqc/fastqc_data.txt" % (c.outdirs[2], sid)
            q22 = "%s/%s_2.PE_fastqc/fastqc_data.txt" % (c.outdirs[2], sid)
            r21, r22 = parse_fastqc(q21), parse_fastqc(q22)
            rc21, rc22 = r21['Total Sequences'], r22['Total Sequences']
            assert rc21 == rc22, "%s: read1 %d != read2 %d" % \
                    (sid, rc21, rc22)
            q31 = "%s/%s_1.SE_fastqc/fastqc_data.txt" % (c.outdirs[2], sid)
            q32 = "%s/%s_2.SE_fastqc/fastqc_data.txt" % (c.outdirs[2], sid)
            r31, r32 = parse_fastqc(q31), parse_fastqc(q32)
            rc31, rc32 = r31['Total Sequences'], r32['Total Sequences']
            t['ReadPairCount'][i] = rc11
            t["TrimmedReadFile1Paired"][i] = f1p
            t["TrimmedReadFile1Unpaired"][i] = f1s
            t["TrimmedReadFile2Paired"][i] = f2p
            t["TrimmedReadFile2Unpaired"][i] = f2s
            t["RetainedReadPairCount"][i] = rc21
            t["RetainedRead1Count"][i] = rc31
            t["RetainedRead2Count"][i] = rc32
    else:
        for i in range(nrow):
            sid, f1 = t['sid'][i], t['Readfile'][i]
            p1 = "%s/%s.fq.gz" % (c.outdirs[1], sid)
            assert op.isfile(p1), "%s not exist" % p1
            pre1 = op.basename(f1).split(".")[0]
            q1 = "%s/%s_fastqc/fastqc_data.txt" % (c.outdirs[0], pre1)
            q2 = "%s/%s_fastqc/fastqc_data.txt" % (c.outdirs[2], sid)
            r1, r2 = parse_fastqc(q1), parse_fastqc(q2)
            rc1, rc2 = r1['Total Sequences'], r2['Total Sequences']
            t['ReadCount'][i] = rc1
            t["TrimmedReadFile"][i] = p1
            t['TrimmedReadCount'][i] = rc2
    t.write(c.olist, format='ascii.tab', overwrite=True)

if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fastq pipelines'
    )
    parser.add_argument('--config', "--cfg", default = "config.ini", help = 'config file')
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("trim", help = "Trim and QC fastq files")
    sp1.add_argument("--check", action = 'store_true', help = "run script in check mode")
    sp1.set_defaults(func = fq_trim)
 
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


#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import numpy as np
import sys
import time
import argparse
import itertools
import luigi

def run_cmds_ssh(cmds):
    username = 'zhou379'
    hostname = 'mesabi.msi.umn.edu'
    port = 22
    pubkey = os.path.join(os.environ['HOME'], '.ssh', 'id_rsa')
    
    key = paramiko.RSAKey.from_private_key_file(pubkey)
    s = paramiko.SSHClient()
    s.load_system_host_keys()
    s.connect(hostname, port, pkey=key)
    for cmd in cmds:
        stdin, stdout, stderr = s.exec_command(cmd)
        for line in stdout:
            print "... " + line.strip("\n")
    s.close()

class WaitJob(luigi.ExternalTask):
    fc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.fc)
class shortread0Sra(luigi.Task):
    pid = luigi.Parameter()
    name = luigi.Parameter(default="shortread0Sra")
    dirs = luigi.Parameter(default="/scratch.global/zhoux379/shortread")
    def run(self):
        pid, name, dirs = self.pid, self.name, self.dirs
        dirw = op.join(dirs, pid)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        if not op.isdir("cps"): os.makedirs("cps")
        assert op.isfile("01.srr.tsv"), "no 01.srr.tsv in %s" % dirw
        ary = np.genfromtxt("01.srr.tsv", names = True, dtype = None, delimiter = "\t")
        cols = list(ary.dtype.names)
        if not op.isdir("03.sra"): os.makedirs("03.sra")
        if not op.isdir("05.reads"): os.makedirs("05.reads")
        d03 = "03.sra" #op.abspath("03.sra")
        d05 = "05.reads" #op.abspath("05.reads")
        fhc = open("04.fastqdump.sh", "w")
        os.chdir("03.sra")
        fnames = []
        for row in ary:
            rid = row[0]
            print rid
            cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra" % (rid[0:3], rid[0:6], rid, rid)
            #os.system(cmd)
            f_sra = "%s.sra" % rid
            assert op.isfile(f_sra), "%s is not there" % f_sra
            fnames.append(f_sra)
            print >>fhc, "fastq-dump --gzip --split-files -outdir %s %s/%s" % (d05, d03, f_sra)
        fname_str = " ".join(fnames)
        fhc.close()
        os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
        jname = name + "Job"
        assert op.isfile(name), "no %s in pbs" % name
        os.system("qsub %s -v JOB=%s,DIR=%s" % (name, jname, dirw))
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        dirc = "%s/%s/cps" % (self.dirs, self.pid)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class shortread2Trim(luigi.Task):
    pid = luigi.Parameter()
    name = luigi.Parameter(default="shortread2Trim")
    dirs = luigi.Parameter(default="/scratch.global/zhoux379/shortread")
    trimm = luigi.Parameter(default="%s/git/trimmomatic" % os.environ["src"])
    def requires(self):
        dirw = op.join(self.dirs, self.pid)
        dirc = "%s/cps" % dirw
        f11 = "%s/11.tsv" % dirw
        if op.isfile(f11) and os.stat(f11).st_size > 0:
            return None
        else:
            return [shortread0Sra(self.pid), WaitJob("%s/shortread0SraJob" % dirc)]
    def run(self):
        pid, name, dirs = self.pid, self.name, self.dirs
        dirw = op.join(dirs, pid)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        if not op.isdir("cps"): os.makedirs("cps")
        assert op.isfile("10.adapter.fas"), "no 10.adapter.fas in %s" % dirw
        assert op.isfile("11.tsv"), "no 11.tsv in %s" % dirw
        ary = np.genfromtxt("11.tsv", names = True, dtype = None, delimiter = "\t")
        f12, f14, f16 = "12.fastqc.sh", "14.trim.sh", "16.fastqc.sh"
        d13, d15, d17 = "13.fastqc", "15.trim", "17.fastqc"
        fho1, fho2, fho3 = open(f12, "w"), open(f14, "w"), open(f16, "w")
        trimm = self.trimm
        assert op.isdir(trimm), "%s is not there" % trimm
        for diro in [d13, d15, d17]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        for row in ary:
            rid, pe, dirf, read1, read2 = list(row)[0:5]
            if pe == 0:
                f1 = "%s/%s" % (dirf, read1)
                assert op.isfile(f1), "%s not there" % f1
                print >>fho1, "fastqc -o %s --noextract -f fastq %s" % (d13, f1)
                fo = "%s/%s" % ("15.trim", read1)
                print >>fho2, "java -jar %s/trimmomatic.jar SE %s %s ILLUMINACLIP:10.adapter.fas:2:30:10:8:no LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % (trimm, f1, fo)
                print >>fho3, "fastqc -o %s --noextract -f fastq %s" % (d17, fo)
            else:
                f1 = "%s/%s" % (dirf, read1)
                f2 = "%s/%s" % (dirf, read2)
                assert op.isfile(f1), "%s not there" % f1
                assert op.isfile(f2), "%s not there" % f2
                pre1, aff1 = read1.split(".")
                pre2, aff2 = read2.split(".")
                f11 = "15.trim/%s.PE.%s" % (pre1, aff1)
                f12 = "15.trim/%s.SE.%s" % (pre1, aff1)
                f21 = "15.trim/%s.PE.%s" % (pre2, aff2)
                f22 = "15.trim/%s.SE.%s" % (pre2, aff2)
                print >>fho1, "fastqc -o %s --noextract -f fastq %s %s" % (d13, f1, f2)
                print >>fho2, "java -jar %s/trimmomatic.jar PE %s %s %s %s %s %s ILLUMINACLIP:10.adapter.fas:2:30:10:8:no LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % (trimm, f1, f2, f11, f12, f21, f22)
                print >>fho3, "fastqc -o %s --noextract -f fastq %s %s %s %s" % (d17, f11, f12, f21, f22)
        os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
        jname = name + "Job"
        assert op.isfile(name), "no %s in pbs" % name
        os.system("qsub %s -v JOB=%s,DIR=%s" % (name, jname, dirw))
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        dirc = "%s/%s/cps" % (self.dirs, self.pid)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))

if __name__ == "__main__":
    luigi.run()

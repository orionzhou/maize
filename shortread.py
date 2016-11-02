#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import re
import numpy as np
import sys
import time
import argparse
import itertools
import pysam
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
def read_config(fcfg):
    cfg = dict()
    for line in open(fcfg, "r"):
        line = line.strip("\n")
        key, val = re.split(' |=|:')
        if key == "paired":
            val = str(val).lower() in ("yes", "true", "t", "1")
        cfg[key] = val
    return cfg
def check_bam(fbam):
    exist_bam = 1
    try:
        bam = pysam.AlignmentFile(fbam, "rb")
    except:
        exist_bam = 0
    return exist_bam

class WaitJob(luigi.ExternalTask):
    fc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.fc)
class shortread1Trim(luigi.Task):
    name = luigi.Parameter(default="shortread1Trim")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    seqlist = luigi.Parameter()
    adapter = luigi.Parameter()
    trimmomatic = luigi.Parameter()
    def run(self):
        name, dirw, seqlist, paired, f_adp, trimm = self.name, self.dirw, self.seqlist, self.paired, self.adapter, self.trimmomatic
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        if not op.isdir("cps"): os.makedirs("cps")
        assert op.isfile(f_adp), "%s not exist" % f_adp 
        assert op.isfile(seqlist), "%s not exist" % seqlist
        ary = np.genfromtxt(seqlist, names = True, dtype = None, delimiter = "\t")
        f12, f14, f16 = "12.fastqc.sh", "14.trim.sh", "16.fastqc.sh"
        d13, d15, d17 = "13.fastqc", "15.trim", "17.fastqc"
        fho1, fho2, fho3 = open(f12, "w"), open(f14, "w"), open(f16, "w")
        assert op.isdir(trimm), "%s is not there" % trimm
        for diro in [d13, d15, d17]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        for row in ary:
            if paired:
                rid, dirf, read1, read2 = list(row)[0:4]
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
                print >>fho2, "java -Xmx2500M -jar %s/trimmomatic.jar PE -threads 4 %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10:8:no LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % (trimm, f1, f2, f11, f12, f21, f22, f_adp)
                print >>fho3, "fastqc -o %s --noextract -f fastq %s %s %s %s" % (d17, f11, f12, f21, f22)
            else:
                rid, dirf, read1 = list(row)[0:3]
                f1 = "%s/%s" % (dirf, read1)
                assert op.isfile(f1), "%s not there" % f1
                print >>fho1, "fastqc -o %s --noextract -f fastq %s" % (d13, f1)
                fo = "%s/%s" % ("15.trim", read1)
                print >>fho2, "java -Xmx2500M -jar %s/trimmomatic.jar SE -threads 4 %s %s ILLUMINACLIP:%s:2:30:10:8:no LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % (trimm, f1, fo, f_adp)
                print >>fho3, "fastqc -o %s --noextract -f fastq %s" % (d17, fo)
        os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
        jname = name + "Job"
        assert op.isfile(name), "no %s in pbs" % name
        os.system("qsub %s -v JOB=%s,DIR=%s" % (name, jname, dirw))
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))
class shortread2Check(luigi.Task):
    name = luigi.Parameter(default="shortread2Check")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    seqlist = luigi.Parameter()
    def requires(self):
        return [shortread1Trim(), WaitJob("%s/cps/shortread1TrimJob" % self.dirw)]
    def run(self):
        name, dirw, paired, seqlist = self.name, self.dirw, self.paired, self.seqlist
        os.chdir(dirw)
        assert op.isfile(seqlist), "%s not exist" % seqlist
        ary = np.genfromtxt(seqlist, names = True, dtype = None, delimiter = "\t")
        cols = list(ary.dtype.names)
        d15, f18 = "15.trim", "18.tsv"
        d15 = op.abspath(d15)
        fho = open(f18, "w")
        if paired:
            print >>fho, "\t".join([cols[0], "dirf", "read1p", "read1s", "read2p", "read2s"]+cols[4:len(cols)])
        else:
            print >>fho, "\t".join([cols[0], "dirf", "read1"]+cols[3:len(cols)])
        for row in ary:
            row = list(row)
            if paired:
                rid, dirf, read1, read2 = row[0:4]
                pre1, aff1 = read1.split(".")
                pre2, aff2 = read2.split(".")
                f1p = "%s.PE.%s" % (pre1, aff1)
                f1s = "%s.SE.%s" % (pre1, aff1)
                f2p = "%s.PE.%s" % (pre2, aff2)
                f2s = "%s.SE.%s" % (pre2, aff2)
                p1p, p1s, p2p, p2s = ["%s/%s" % (d15, x) for x in [f1p, f1s, f2p, f2s]]
                assert op.isfile(p1p), "%s not exist" % p1p
                assert op.isfile(p1s), "%s not exist" % p1s
                assert op.isfile(p2p), "%s not exist" % p2p
                assert op.isfile(p2s), "%s not exist" % p2s
                print >>fho, "\t".join([rid, d15, f1p, f1s, f2p, f2s] + row[4:len(row)])
            else:
                rid, dirf, read1 = row[0:3]
                f1 = read1
                p1 = "%s/%s" % (d15, f1)
                assert op.isfile(p1), "%s not exist" % p1
                print >>fho, "\t".join([rid, d15, f1]+row[3:len(row)])
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))

class shortread3Tophat(luigi.Task):
    name = luigi.Parameter(default="shortread3Tophat")
    species = luigi.Parameter()
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    seqlist = luigi.Parameter(default="18.tsv")
    db = luigi.Parameter()
    samstat = luigi.Parameter()
    def run(self):
        name, species, dirw, paired, seqlist, db, samstat = self.name, self.species, self.dirw, self.paired, self.seqlist, self.db, self.samstat
        os.chdir(dirw)
        assert op.isfile(seqlist), "%s not exist" % seqlist
        assert op.isfile(samstat), "%s not exist" % samstat
        ary = np.genfromtxt(seqlist, names = True, dtype = None, delimiter = "\t")
        f211, f212, f213 = "21.1.tophat.sh", "21.2.index.sh", "21.3.stat.sh"
        d22 = "22.tophat"
        fho1, fho2, fho3 = open(f211, "w"), open(f212, "w"), open(f213, "w")
        for diro in [d22]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        min_intron_len = 60000
        if species == 'rice': min_intron_len = 20000
        for row in ary:
            row = list(row)
            fbam = "%s/%s/accepted_hits.bam" % (d22, row[0])
            exist_fbam = check_bam(fbam)
            if paired:
                rid, dirf, read1p, read1s, read2p, read2s = row[0:6]
                f1 = "%s/%s" % (dirf, read1p)
                f2 = "%s/%s" % (dirf, read2p)
                if not exist_fbam:
                    print >>fho1, "tophat2 --num-threads 24 --max-multihits 20 --min-intron-length 5 --max-intron-length %d -o %s/%s %s %s %s" % (min_intron_len, d22, rid, db, f1, f2)
                print >>fho2, "samtools index %s/%s/accepted_hits.bam" % (d22, rid)
                print >>fho3, "%s %s/%s/accepted_hits.bam" % (samstat, d22, rid)
            else:
                rid, dirf, read1 = row[0:3]
                f1 = "%s/%s" % (dirf, read1)
                if not exist_fbam:
                    print >>fho1, "tophat2 --num-threads 24 --max-multihits 20 --min-intron-length 5 --max-intron-length %d -o %s/%s %s %s" % (min_intron_len, d22, rid, db, f1)
                print >>fho2, "samtools index %s/%s/accepted_hits.bam" % (d22, rid)
                print >>fho3, "%s %s/%s/accepted_hits.bam" % (samstat, d22, rid)
        os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
        jname = name + "Job"
        assert op.isfile(name), "no %s in pbs" % name
        os.system("qsub %s -v JOB=%s,DIR=%s" % (name, jname, dirw))
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))
class shortread4Htseq(luigi.Task):
    name = luigi.Parameter(default="shortread4Htseq")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    stranded = luigi.BoolParameter()
    bamlist = luigi.Parameter(default="23.tsv")
    annotation = luigi.Parameter()
    def run(self):
        name, dirw, paired, stranded, bamlist, annotation = self.name, self.dirw, self.paired, self.stranded, self.bamlist, self.annotation
        os.chdir(dirw)
        assert op.isfile(bamlist), "%s not exist" % bamlist
        assert op.isfile(annotation), "%s not exist" % annotation
        ary = np.genfromtxt(bamlist, names = True, dtype = None, delimiter = "\t")
        f31 = "31.htseq.sh"
        d32 = "32.htseq"
        fho1 = open(f31, "w")
        for diro in [d32]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        for row in ary:
            row = list(row)
            rid, fbam = row[0:2]
            print >>fho1, "htseq-count -s no -t gene -i ID -m union -a 20 -f bam %s %s > %s/%s.txt" % (fbam, annotation, d32, rid)
        os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
        jname = name + "Job"
        assert op.isfile(name), "no %s in pbs" % name
        #os.system("qsub %s -v JOB=%s,DIR=%s" % (name, jname, dirw))
        #os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))

if __name__ == "__main__":
    luigi.run()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import time
import argparse
import paramiko
import luigi

def run_cmds_ssh(cmds):
    username = 'zhoup'
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
class GenomeFas(luigi.Task):
    org = luigi.Parameter()
    name = luigi.Parameter(default="GenomeFas")
    species = luigi.Parameter(default="maize")
    dirg = luigi.Parameter(default="/home/springer/zhoux379/Data/genome")
    dird = luigi.Parameter(default="/home/springer/zhoux379/Data/db")
    ref = luigi.IntParameter(default=0)
    def run(self):
        org, name, dirg, dird = self.org, self.name, self.dirg, self.dird
        dirw = op.join(dirg, org)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        if not op.isdir("cps"): os.makedirs("cps")

        for fname in ["raw.fix.fas.index", "11_genome.fas.index"]:
            if op.isfile(fname):
                os.remove(fname)
        if op.islink("11_genome.fas"): os.unlink("11_genome.fas")
        
        assert op.isfile("raw.fas"), "no raw.fas in %s" % dirw
        #os.system("seq.check.pl -i raw.fas -o raw.fix.fas")
        if self.ref == 1:
            print self.ref
            #os.system("cp raw.fix.fas 11_genome.fas")
        else:
            os.system("seq.rename.pl -i raw.fix.fas -p scf -o 11_genome.fas")
            if op.isfile("ctg.raw.fas"):
                os.system("seq.check.pl -i ctg.raw.fas -o ctg.fas")
        fg = "%s/11_genome.fas" % dirw
        #os.system("seqlen.py 11_genome.fas 15.sizes")
        #os.system("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print $1, 0, $2}' 15.sizes > 15.bed")
        #os.system("seqgap.pl -i 11_genome.fas -o 16.gap.bed -m 10")
        os.system("bedToBigBed 16.gap.bed 15.sizes 16.gap.bb")
        os.system("bgzip -c 16.gap.bed > 16.gap.bed.gz")
        os.system("tabix -p bed 16.gap.bed.gz")
        
        if op.isdir("12.rm"): os.system("rm -rf 12.rm")
        os.makedirs("12.rm")
        jname = name + "Job"
        os.chdir("%s/pbs" % op.dirname(op.realpath(__file__)))
        #os.system("qsub %s -v DIR=%s,SPE=%s,JOB=%s" % (name, dirw, self.species, jname))
        
        if not op.isdir(dird): os.makedirs(dird)
        dird1 = op.join(dird, "blat")
        if not op.isdir(dird1): os.makedirs(dird1)
        os.chdir(dird1)
        #os.system("module load blat")
        os.system("faToTwoBit %s %s.2bit" % (fg, org))
        os.system("blat %s.2bit tmp.fas tmp.out -makeOoc=%s.2bit.tile11.ooc" % (org, org))
        if op.isfile("tmp.out"): os.remove("tmp.out")

        dird2 = op.join(dird, 'bowtie2')
        if not op.isdir(dird2): os.makedirs(dird2)
        os.chdir(dird2)
        fd = "%s.fas" % org
        if op.isfile(fd): os.system("rm %s" % fd)
        fd = "%s.fa" % org
        os.system("ln -sf %s %s" % (fg, fd))
        #os.system("module load bowtie2")
        os.system("bowtie2-build %s %s" % (fd, org))
        
        dirc = "%s/%s/cps" % (self.dirg, self.org)
        os.system("touch %s/%s" % (dirc, name))
    def output(self):
        dirc = "%s/%s/cps" % (self.dirg, self.org)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class RepeatMaskerParse(luigi.Task):
    org = luigi.Parameter()
    name = luigi.Parameter(default="RepeatMaskerParse")
    def requires(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return [RepeatMasker(self.org), WaitJob("%s/RepeatMaskerJob" % dirc)]
    def run(self):
        org, name, = self.org, self.name
        dirw = op.join(os.environ['genome'], org)
        os.chdir(dirw)
        os.system("parse.rm.pl -i 12.rm/11_genome.fas.out -o 12.rm.tbl")
        os.system("awk 'BEGIN{OFS=\"\\t\"} {print $1, $2-1, $3, $9 \" | \" ($5)}' 12.rm.tbl > 12.rm.raw.bed")
        os.system("sortBed -i 12.rm.raw.bed > 12.rm.bed")
        os.system("bedToBigBed -tab -type=bed4 12.rm.bed 15.sizes 12.rm.bb")
        os.system("rm 12.rm.raw.bed")
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))

if __name__ == "__main__":
    luigi.run()

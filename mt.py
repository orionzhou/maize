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
    def run(self):
        org = self.org
        dirw = op.join(os.environ['genome'], org)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        if not op.isdir("cps"): os.makedirs("cps")

        for fname in ["raw.fix.fas.index", "11_genome.fas.index"]:
            if op.isfile(fname):
                os.remove(fname)
        if op.islink("11_genome.fas"): os.unlink("11_genome.fas")
        
        assert op.isfile("raw.fas"), "no raw.fas in %s" % dirw
        os.system("seq.check.pl -i raw.fas -o raw.fix.fas")
        os.system("seq.rename.pl -i raw.fix.fas -p scf -o 11_genome.fas")
        if op.isfile("ctg.raw.fas"):
            os.system("seq.check.pl -i ctg.raw.fas -o ctg.fas")
        os.system("seqlen.py 11_genome.fas 15.sizes")
        os.system("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print $1, 0, $2}' 15.sizes > 15.bed")
        os.system("seqgap.pl -i 11_genome.fas -o 16.gap.bed -m 10")
        os.system("bedToBigBed 16.gap.bed 15.sizes 16.gap.bb")
        os.system("bgzip -c 16.gap.bed > 16.gap.bed.gz")
        os.system("tabix -p bed 16.gap.bed.gz")
        os.system("touch cps/GenomeFas")
    def output(self):
        return luigi.LocalTarget(op.join(os.environ['genome'], self.org, "cps", "GenomeFas"))
class GenomeDb(luigi.Task):
    org = luigi.Parameter()
    def requires(self):
        return GenomeFas(self.org)
    def run(self):
        org = self.org
        dirw = op.join(os.environ['genome'], org)
        fg = op.join(dirw, "11_genome.fas")
        assert op.isfile(fg), "no 11_genome.fas found in %s" % dirw
        
        dird = op.join(os.environ['data'], 'db', 'blat')
        os.chdir(dird)
        os.system("faToTwoBit %s %s.2bit" % (fg, org))
        os.system("blat %s.2bit tmp.fas tmp.out -makeOoc=%s.2bit.tile11.ooc" % (org, org))
        if op.isfile("tmp.out"): os.remove("tmp.out")

        dird = op.join(os.environ['data'], 'db', 'bowtie2')
        os.chdir(dird)
        fd = "%s.fas" % org
        if op.isfile(fd): os.system("rm %s" % fd)
        fd = "%s.fa" % org
        os.system("ln -sf %s %s" % (fg, fd))
        os.system("bowtie2-build %s %s" % (fd, org))
        os.system("touch %s/%s/cps/GenomeDb" % (os.environ['genome'], org))
    def output(self):
        return luigi.LocalTarget(op.join(os.environ['genome'], self.org, "cps", "GenomeDb"))
class RepeatMasker(luigi.Task):
    org = luigi.Parameter()
    name = luigi.Parameter(default="RepeatMasker")
    def requires(self):
        return GenomeFas(self.org)
    def run(self):
        org, name, = self.org, self.name
        dirw = op.join(os.environ['genome'], org)
        os.chdir(dirw)
        fg = op.join(dirw, "11_genome.fas")
        if op.isdir("12.rm"): os.system("rm -rf 12.rm")
        os.makedirs("12.rm")
        nname = name + "Job"
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        os.chdir("%s/pbs" % os.environ['code'])
        os.system("qsub %s -v ORG=%s" % (name, org))
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
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

class Comp1Pre(luigi.Task):
    qry = luigi.Parameter()
    tgt = luigi.Parameter(default="HM101")
    name = luigi.Parameter(default="Comp1Pre")
    def requires(self):
        return [GenomeDb(self.qry), RepeatMaskerParse(self.qry)]
    def run(self):
        qry, tgt, name = self.qry, self.tgt, self.name
        dird = os.environ['data']
        qry_fas = "%s/genome/%s/11_genome.fas" % (dird, qry)
        tgt_fas = "%s/genome/%s/11_genome.fas" % (dird, tgt)
        qry_2bit = "%s/db/blat/%s.2bit" % (dird, qry)
        tgt_2bit = "%s/db/blat/%s.2bit" % (dird, tgt)
        qry_size = "%s/genome/%s/15.sizes" % (dird, qry)
        tgt_size = "%s/genome/%s/15.sizes" % (dird, tgt)
        qry_size_bed = "%s/genome/%s/15.bed" % (dird, qry)
        tgt_size_bed = "%s/genome/%s/15.bed" % (dird, tgt)
        qry_gap = "%s/genome/%s/16.gap.bed" % (dird, qry)
        tgt_gap = "%s/genome/%s/16.gap.bed" % (dird, tgt)
        dirw = "%s/misc3/%s_%s/23_blat" % (dird, qry, tgt)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        dirc = "%s/misc3/%s_%s/cps" % (dird, qry, tgt)
        if not op.isdir(dirc): os.makedirs(dirc)

        if not op.isdir("01_seq"): os.makedirs("01_seq")
        os.system("breakseq.bygap.pl -i %s -o 00.fas -g 1000" % qry_fas)
        os.system("seq.splitlarge.py 00.fas 00.even.fas")
        os.system("qsub.blat.pl -i 00.even.fas -o 01_seq -n 10 -t %s -g %s" % (tgt, qry))
        
        os.chdir("%s/pbs" % os.environ['code'])
        nname = 'Comp2Blat'
        for i in range(0, 10):
            cmd = "qsub %s -N %s.%d -v JOB=%s.%d,PRE=%s/01_seq/part,SUF=fas,PPN=%s,BAT=%d,DIG=3,QRY=%s,TGT=%s" % (nname, nname, i, nname, i, dirw, 24, i, qry, tgt)
            os.system(cmd)
        ### waiting for Comp2Blat.*
        
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Comp3Process(luigi.Task):
    qry = luigi.Parameter()
    tgt = luigi.Parameter(default="HM101")
    name = luigi.Parameter(default="Comp3Process")
    def requires(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return [Comp1Pre(self.qry, self.tgt)] + [WaitJob("%s/Comp2Blat.%d" % (dirc, i)) for i in range(0,10)]
    def run(self):
        qry, tgt, name = self.qry, self.tgt, self.name
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        nname = name + "Job"
        os.chdir("%s/pbs" % os.environ['code'])
        os.system("qsub %s -v JOB=%s,QRY=%s,TGT=%s" % (name, nname, qry, tgt))
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Comp4Blatnov(luigi.Task):
    qry = luigi.Parameter()
    tgt = luigi.Parameter(default="HM101")
    name = luigi.Parameter(default="Comp4Blatnov")
    def requires(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return [Comp3Process(self.qry, self.tgt), WaitJob("%s/Comp3ProcessJob" % dirc)]
    def run(self):
        qry, tgt, name = self.qry, self.tgt, self.name
        dirw = "%s/misc3/%s_%s/23_blat" % (os.environ['data'], qry, tgt)
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], qry, tgt)
        nname = name + "Job"
        os.chdir("%s/pbs" % os.environ['code'])
        cmd = "qsub %s -v JOB=%s,PRE=%s/24.nov/part,SUF=fas,PPN=%s,BAT=%d,DIG=2,QRY=%s,TGT=%s" % (name, nname, dirw, 24, 0, qry, tgt)
        os.system(cmd)
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Comp5Process(luigi.Task):
    qry = luigi.Parameter()
    tgt = luigi.Parameter(default="HM101")
    name = luigi.Parameter(default="Comp5Process")
    def requires(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return [Comp4Blatnov(self.qry, self.tgt), WaitJob("%s/Comp4BlatnovJob" % dirc)]
    def run(self):
        qry, tgt, name, = self.qry, self.tgt, self.name
        nname = name + "Job"
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        os.chdir("%s/pbs" % os.environ['code'])
        os.system("qsub %s -v JOB=%s,QRY=%s,TGT=%s" % (name, nname, qry, tgt))
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Comp6Novseq(luigi.Task):
    qry = luigi.Parameter()
    tgt = luigi.Parameter(default="HM101")
    name = luigi.Parameter(default="Comp6Novseq")
    def requires(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return [Comp5Process(self.qry, self.tgt), WaitJob("%s/Comp5ProcessJob" % dirc)]
    def run(self):
        qry, tgt, name = self.qry, self.tgt, self.name
        dird = os.environ['data']
        qry_fas = "%s/genome/%s/11_genome.fas" % (dird, qry)
        tgt_fas = "%s/genome/%s/11_genome.fas" % (dird, tgt)
        qry_2bit = "%s/db/blat/%s.2bit" % (dird, qry)
        tgt_2bit = "%s/db/blat/%s.2bit" % (dird, tgt)
        qry_size = "%s/genome/%s/15.sizes" % (dird, qry)
        tgt_size = "%s/genome/%s/15.sizes" % (dird, tgt)
        qry_size_bed = "%s/genome/%s/15.bed" % (dird, qry)
        tgt_size_bed = "%s/genome/%s/15.bed" % (dird, tgt)
        qry_gap = "%s/genome/%s/16.gap.bed" % (dird, qry)
        tgt_gap = "%s/genome/%s/16.gap.bed" % (dird, tgt)
        dirw = "%s/misc3/%s_%s/41_novseq" % (dird, qry, tgt)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        os.system("gax2bed.pl -i ../23_blat/41.5/gax -p tgt -o - | sortBed -i stdin | mergeBed -i stdin > 00.bed")
        os.system("subtractBed -a %s -b %s | subtractBed -a stdin -b 00.bed | bedfilter.pl -l 50 -o 01.bed" % (qry_size_bed, qry_gap))
        os.system("seqret.pl -d %s -b 01.bed -o 01.fas" % qry_fas)
        dirc = "%s/misc3/%s_%s/cps" % (dird, qry, tgt)
        if not op.isdir(dirc): os.makedirs(dirc)
        nname = name + "Job"
        os.chdir("%s/pbs" % os.environ['code'])
        os.system("qsub %s -v JOB=%s,QRY=%s,TGT=%s" % (name, nname, qry, tgt))
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Comp6Novseq2(luigi.Task):
    qry = luigi.Parameter()
    tgt = luigi.Parameter(default="HM101")
    name = luigi.Parameter(default="Comp6Novseq2")
    def requires(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return [Comp6Novseq(self.qry, self.tgt), WaitJob("%s/Comp6NovseqJob" % dirc)]
    def run(self):
        qry, tgt, name, = self.qry, self.tgt, self.name
        dird = os.environ['data']
        qry_fas = "%s/genome/%s/11_genome.fas" % (dird, qry)
        dirw = "%s/misc3/%s_%s/41_novseq" % (os.environ['data'], qry, tgt)
        os.chdir(dirw)
        os.system("blastnr.pl -i 01.fas -o 11.blastnr")
        os.system("comp.novseq.pl -q %s -t %s" % (qry, tgt))
        os.system("awk 'BEGIN{OFS=\"\\t\"} {if($4==\"foreign\") print}' 12.bed > 12.foreign.bed")
        os.system("awk 'BEGIN{OFS=\"\\t\"} {if($4!=\"foreign\") print}' 12.bed > 12.unforeign.bed")
        os.system("ln -sf 12.unforeign.bed 21.bed")
        os.system("seqret.pl -d %s -b 21.bed -o 21.fas" % qry_fas)
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/misc3/%s_%s/cps" % (os.environ['data'], self.qry, self.tgt)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))

class Anno1Pre(luigi.Task):
    org = luigi.Parameter()
    name = luigi.Parameter(default="Anno1Pre")
    def requires(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return Comp6Novseq2(self.org)
    def run(self):
        org, name, = self.org, self.name
        dirw = op.join(os.environ['genome'], org, 'augustus')
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        assert op.isfile("21.gff"), "no 21.gff in %s" % dirw
        assert op.isfile("../raw.fix.fas.map"), "no ../raw.fix.fas.map in %s" % dirw
        os.system("gff.ncgr.py --map ../raw.fix.fas.map 21.gff 22.gff")
        os.system("gff2gtb.pl -i 22.gff -o 22.gtb")
        os.system("gtb.rmutr.pl -i 22.gtb -o 23.gtb")
        os.system("gtb.dedup.pl -i 23.gtb -o 25.dedup.gtb")
        os.system("gtb2gtbx.pl -i 25.dedup.gtb -d ../11_genome.fas -o 31.gtbx")
        os.system("cut -f1-18 31.gtbx > 31.gtb")
        os.system("gtbx2fas.pl -i 31.gtbx -o 31.fas")
        os.system("gtb2gff.pl -i 31.gtb -o 31.gff")
        os.chdir("%s/pbs" % os.environ['code'])
        nname = name + "Job"
        os.system("qsub %s -v JOB=%s,ORG=%s" % (name, nname, org))
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Anno2Pfam(luigi.Task):
    org = luigi.Parameter()
    name = luigi.Parameter(default="Anno2Pfam")
    def requires(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return [Anno1Pre(self.org), WaitJob("%s/Anno1PreJob" % dirc)]
    def run(self):
        org, name, = self.org, self.name
        dirw = op.join(os.environ['genome'], org, 'augustus')
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        assert op.isfile("34.1.txt"), "no 34.1.txt in %s" % dirw
        f_pfam = "%s/db/pfam/Pfam-A.hmm" % os.environ['data']
        os.system("hmmc2htb.pl -i 34.1.txt -o 34.2.htb -m %s -s 31.fas" % f_pfam)
        os.system("htb.qtile.pl -i 34.2.htb -o 34.3.htb")
        os.system("htb.filter.pl -i 34.3.htb -l 10 -e 0.01 -o 34.4.htb")
        os.system("cut -f2-4,6,7-9,11-13 34.4.htb > 34.tsv")
        os.system("gtb.addpfam.pl -i 31.gtb -p 34.tsv -o 41.gtb")
        os.system("gtb2bed.s.pl -i 41.gtb -o 41.bed")
        assert op.isfile("../12.rm.bed"), "no ../12.rm.bed in %s" % dirw
        os.system("intersectBed -wao -a 41.bed -b ../12.rm.bed > 42.ovlp.bed")
        os.system("gtb.addrm.pl -i 41.gtb -v 42.ovlp.bed -o 43.gtb")
        os.system("mt.nbs.pl -g %s" % org)
        os.system("mt.rlk.pl -g %s" % org)
        os.chdir("%s/pbs" % os.environ['code'])
        nname = name + "Job"
        os.system("qsub %s -v JOB=%s,ORG=%s" % (name, nname, org))
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))
class Anno3Post(luigi.Task):
    org = luigi.Parameter()
    name = luigi.Parameter(default="Anno3Post")
    def requires(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return [Anno2Pfam(self.org), WaitJob("%s/Anno2PfamJob" % dirc)]
    def run(self):
        org, name, = self.org, self.name
        dirw = op.join(os.environ['genome'], org)
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        assert op.isfile("augustus/43.gtb"), "no augustus/43.gtb in %s" % dirw
        os.system("ln -sf augustus/43.gtb 41.gtb")
        assert op.isfile("42.nbs/11.gtb"), "no 42.nbs/11.gtb in %s" % dirw
        os.system("ln -sf 42.nbs/11.gtb 42.nbs.gtb")
        assert op.isfile("43.crp/61_final.gtb"), "no 43.crp/61_final.gtb in %s" % dirw
        os.system("gtb.crp.pl -i 43.crp/61_final.gtb -o 43.crp.gtb")
        assert op.isfile("44.rlk/11.gtb"), "no 44.rlk/11.gtb in %s" % dirw
        os.system("ln -sf 44.rlk/11.gtb 44.rlk.gtb")
        os.system("gtb.merge.pl -a 41.gtb -b 42.nbs.gtb -o 49.1.gtb")
        os.system("gtb.merge.pl -a 49.1.gtb -b 43.crp.gtb -o 49.2.gtb")
        os.system("gtb.merge.pl -a 49.2.gtb -b 44.rlk.gtb -o 49.gtb")
        os.system("gtb.dedup.pl -i 49.gtb -o 50.1.dedup.gtb")
        os.system("gtb.pickalt.pl -i 50.1.dedup.gtb -o 50.2.pickalt.gtb")
        f_ctm = "%s/%s_HM101/41_novseq/15.foreign.scf.txt"
        if org == "HM101" or not op.isfile(f_ctm):
            print "%s not there: skip contaminant removal" % f_ctm
            os.system("gtb.fill.len.pl -i 50.2.pickalt.gtb | gtb.filter.pl -l 30 -o 51.gtb")
        else:
            os.system("gtb.fill.len.pl -i 50.2.pickalt.gtb | gtb.filter.pl -l 30 -c %s -o 51.gtb" % f_ctm)
        os.system("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if(NR==1 || tolower($16) != \"te\") print}' 51.gtb > 55_noTE.gtb")
        os.system("gtb2gff.pl -i 51.gtb -o 51.gff");
        os.system("gtb.idx.pl -i 51.gtb -s 15.sizes");
        os.system("gtb2tbl.pl -i 51.gtb -o 51.tbl");
        os.system("gtb2fas.pl -i 51.gtb -d 11_genome.fas -o 51.fas");
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        os.system("touch %s/%s" % (dirc, self.name))
    def output(self):
        dirc = "%s/%s/cps" % (os.environ['genome'], self.org)
        return luigi.LocalTarget("%s/%s" % (dirc, self.name))

if __name__ == "__main__":
    luigi.run()

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
import luigi
from fadapa import Fadapa
from string import Template
import pysam
from colorama import init, Fore, Back, Style

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
            print("... " + line)
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
def parse_fastqc(fqc):
    assert op.isfile(fqc), "%s not exist" % fqc
    qc = Fadapa(fqc)
    r = dict()
    for ary in qc.clean_data('Basic Statistics'):
        r[ary[0]] = ary[1]
    return r
def parse_samtools_stat(fi):
    assert op.isfile(fi), "%s not exist" % fi
    r = dict()
    for line in open(fi, "r"):
        line = line.strip("\n")
        row = line.split("\t")
        if not row[0] == 'SN':
            continue
        r[row[1].strip(":")] = row[2]
    return r
def parse_tophat_alignsum(fi):
    assert op.isfile(fi), "%s not exist" % fi
    r = dict()
    for line in open(fi, "r"):
        row = line.split(":")
        if len(row) < 2:
            continue
        key = row[0].strip(" \n")
        val = row[1].strip(" \n").split(" ")[0]
        r[key] = val
    return r

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
    pbs_template = luigi.Parameter()
    pbs_queue = luigi.Parameter()
    pbs_walltime = luigi.Parameter()
    pbs_ppn = luigi.Parameter()
    pbs_email = luigi.Parameter()
    def run(self):
        name, dirw, seqlist, paired, f_adp, trimm = self.name, self.dirw, self.seqlist, self.paired, self.adapter, self.trimmomatic
        if not op.isdir(dirw): os.makedirs(dirw)
        os.chdir(dirw)
        if not op.isdir("cps"): os.makedirs("cps")
        assert op.isfile(f_adp), "%s not exist" % f_adp 
        assert op.isfile(seqlist), "%s not exist" % seqlist
        ary = np.genfromtxt(seqlist, names = True, dtype = object, delimiter = "\t")
        fo1, fo2, fo3 = "12.1.fastqc.sh", "12.2.trim.sh", "12.3.fastqc.sh"
        d13, d14, d15 = "13.fastqc", "14.trim", "15.fastqc"
        fho1, fho2, fho3 = open(fo1, "w"), open(fo2, "w"), open(fo3, "w")
        assert op.isfile(trimm), "%s is not there" % trimm
        for diro in [d13, d14, d15]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        for row in ary:
            row = [str(x, 'utf-8') for x in list(row)]
            sid = row[0]
            if paired:
                f1, f2 = row[5:7]
                assert op.isfile(f1), "%s not there" % f1
                assert op.isfile(f2), "%s not there" % f2
                f11 = "%s/%s_1.PE.fastq.gz" % (d14, sid)
                f12 = "%s/%s_1.SE.fastq.gz" % (d14, sid)
                f21 = "%s/%s_2.PE.fastq.gz" % (d14, sid)
                f22 = "%s/%s_2.SE.fastq.gz" % (d14, sid)
                fho1.write("fastqc -o %s --extract -f fastq %s %s\n" % (d13, f1, f2))
                fho2.write("java -Xmx2500M -jar %s PE -threads 4 %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10:8:no LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n" % (trimm, f1, f2, f11, f12, f21, f22, f_adp))
                fho3.write("fastqc -o %s --extract -f fastq %s %s %s %s\n" % (d15, f11, f12, f21, f22))
            else:
                f1 = row[5]
                assert op.isfile(f1), "%s not there" % f1
                fho1.write("fastqc -o %s --noextract -f fastq %s\n" % (d13, f1))
                fo = "%s/%s.fastq.gz" % (d14, sid)
                fob = "%s/%s_1.fastq.gz" % (d14, sid)
                if op.isfile(fob):
                    os.system("mv %s %s" % (fob, fo))
                fho2.write("java -Xmx2500M -jar %s SE -threads 4 %s %s ILLUMINACLIP:%s:2:30:10:8:no LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n" % (trimm, f1, fo, f_adp))
                fho3.write("fastqc -o %s --noextract -f fastq %s\n" % (d15, fo))
        
        cmds = []
        cmds.append("module load fastqc/0.11.2")
        cmds.append("cd %s" % dirw)
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo1))
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo2))
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo3))
        cmds.append("touch cps/%s" % (name + "Job"))
        cmd = "\n".join(cmds)
        
        temdict = {
                "queue": self.pbs_queue,
                "walltime": self.pbs_walltime,
                "ppn": self.pbs_ppn,
                "email": self.pbs_email,
                "cmds": cmd
        }
        fo = "12.pbs"
        fho = open(fo, "w")
        assert op.isfile(self.pbs_template), "cannot read template: %s" % self.pbs_template
        fht = open(self.pbs_template, "r")
        src = Template(fht.read())
        fho.write(src.substitute(temdict))
        
        init()
        print(Fore.GREEN)
        print("A job script has been generated: %s" % fo)
        print("Please check, make necessary changes , then type:")
        print(Fore.RED + "qsub %s" % fo)
        print(Style.RESET_ALL)
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))
class shortread2Check(luigi.Task):
    name = luigi.Parameter(default="shortread2Check")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    ilist = luigi.Parameter(default="00.1.read.tsv")
    olist = luigi.Parameter(default="00.2.trim.tsv")
    def requires(self):
        return [shortread1Trim(), WaitJob("%s/cps/shortread1TrimJob" % self.dirw)]
    def run(self):
        name, dirw, paired, ilist, outseqlist = self.name, self.dirw, self.paired, self.seqlist, self.olist
        os.chdir(dirw)
        assert op.isfile(ilist), "%s not exist" % ilist
        ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
        cols = list(ary.dtype.names)
        do1, do2, do3 = "13.fastqc", "14.trim", "15.fastqc"
        fho = open(olist, "w")
        if paired:
            fho.write("\t".join(cols + ["ReadPairCount", "TrimmedReadFile1Paired", "TrimmedReadFile1Unpaired", "TrimmedReadFile2Paired", "TrimmedReadFile2Unpaired", "RetainedReadPairCount", "RetainedRead1Count", "RetainedRead2Count"]) + "\n")
        else:
            fho.write("\t".join(cols + ["ReadCount", "TrimmedReadFile", "TrimmedReadCount"]) + "\n")
        for row in ary:
            row = [str(x, 'utf-8') for x in list(row)]
            sid = row[0]
            if paired:
                f1, f2 = row[5:7]
                f1p = "%s_1.PE.fastq.gz" % (sid)
                f1s = "%s_1.SE.fastq.gz" % (sid)
                f2p = "%s_2.PE.fastq.gz" % (sid)
                f2s = "%s_2.SE.fastq.gz" % (sid)
                p1p, p1s, p2p, p2s = [op.abspath("%s/%s" % (do2, x)) for x in [f1p, f1s, f2p, f2s]]
                assert op.isfile(p1p), "%s not exist" % p1p
                assert op.isfile(p1s), "%s not exist" % p1s
                assert op.isfile(p2p), "%s not exist" % p2p
                assert op.isfile(p2s), "%s not exist" % p2s
                pre1 = op.basename(f1).rstrip(".fastq.gz")
                pre2 = op.basename(f1).rstrip(".fastq.gz")
                q11 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre1)
                q12 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre2)
                r11, r12 = parse_fastqc(q11), parse_fastqc(q12)
                rc11, rc12 = r11['Total Sequences'], r12['Total Sequences']
                assert rc11 == rc12, "%s: read1 %d != read2 %d" % (sid, rc11, rc12)
                q21 = "%s/%s_1.PE_fastqc/fastqc_data.txt" % (do3, sid)
                q22 = "%s/%s_2.PE_fastqc/fastqc_data.txt" % (do3, sid)
                r21, r22 = parse_fastqc(q21), parse_fastqc(q22)
                rc21, rc22 = r21['Total Sequences'], r22['Total Sequences']
                assert rc21 == rc22, "%s: read1 %d != read2 %d" % (sid, rc21, rc22)
                q31 = "%s/%s_1.SE_fastqc/fastqc_data.txt" % (do3, sid)
                q32 = "%s/%s_2.SE_fastqc/fastqc_data.txt" % (do3, sid)
                r31, r32 = parse_fastqc(q31), parse_fastqc(q32)
                rc31, rc32 = r31['Total Sequences'], r32['Total Sequences']
                fho.write("\t".join(row + [rc11, p1p, p1s, p2p, p2s, rc21, rc31, rc32]) + "\n")
            else:
                f1 = row[5]
                p1 = op.abspath("%s/%s.fastq.gz" % (do2, sid))
                assert op.isfile(p1), "%s not exist" % p1
                pre1 = op.basename(f1).rstrip(".fastq.gz")
                q1 = "%s/%s_fastqc/fastqc_data.txt" % (do1, pre1)
                q2 = "%s/%s_fastqc/fastqc_data.txt" % (do3, sid)
                r1, r2 = parse_fastqc(q1), parse_fastqc(q2)
                rc1, rc2 = r1['Total Sequences'], r2['Total Sequences']
                fho.write("\t".join(row + [rc1, p1, rc2]) + "\n")
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))

class shortread3Tophat(luigi.Task):
    name = luigi.Parameter(default="shortread3Tophat")
    species = luigi.Parameter()
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    ilist = luigi.Parameter(default="00.2.trim.tsv")
    db_bowtie2 = luigi.Parameter()
    tophat2 = luigi.Parameter()
    samstat = luigi.Parameter()
    pbs_template = luigi.Parameter()
    pbs_queue = luigi.Parameter()
    pbs_walltime = luigi.Parameter()
    pbs_ppn = luigi.Parameter()
    pbs_email = luigi.Parameter()
    def requires(self):
        return shortread2Check()
    def run(self):
        name, species, dirw, paired, ilist, db_bowtie2 = self.name, self.species, self.dirw, self.paired, self.ilist, self.db_bowtie2
        os.chdir(dirw)
        assert op.isfile(ilist), "%s not exist" % ilist
        assert op.isfile(self.samstat), "%s not exist" % self.samstat
        ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
        fo1, fo2, fo3, fo4 = "21.1.tophat.sh", "21.2.index.sh", "21.3.samtools.sh", "21.4.samstat.sh"
        d22 = "22.tophat"
        fho1, fho2, fho3, fho4 = open(fo1, "w"), open(fo2, "w"), open(fo3, "w"), open(fo4, "w")
        for diro in [d22]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        min_intron_len = 60000
        if species == 'rice': min_intron_len = 20000
        for row in ary:
            row = [str(x, 'utf-8') for x in list(row)]
            sid = row[0]
            fbam = "%s/%s/accepted_hits.bam" % (d22, row[0])
            exist_fbam = check_bam(fbam)
            if paired:
                f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
                if not exist_fbam:
                    fho1.write("tophat2 --num-threads 24 --max-multihits 20 --min-intron-length 5 --max-intron-length %d -o %s/%s %s %s %s\n" % (min_intron_len, d22, sid, db_bowtie2, f1p, f2p))
            else:
                fr, rc, ft, rrc = row[5:9]
                if not exist_fbam:
                    fho1.write("%s --num-threads 24 --max-multihits 20 --min-intron-length 5 --max-intron-length %d -o %s/%s %s %s\n" % (self.tophat2, min_intron_len, d22, sid, db_bowtie2, ft))
            fho2.write("samtools index %s/%s/accepted_hits.bam\n" % (d22, sid))
            fho3.write("samtools stats %s/%s/accepted_hits.bam > %s/%s/samtools.stat\n" % (d22, sid, d22, sid))
            fho4.write("%s %s/%s/accepted_hits.bam\n" % (self.samstat, d22, sid))

        cmds = []
        cmds.append("module load bowtie2/2.2.4")
        cmds.append("module load samtools")
        cmds.append("module load parallel")
        cmds.append("cd %s" % dirw)
        cmds.append("bash %s" % fo1)
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo2))
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo3))
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo4))
        cmds.append("touch cps/%s" % (name + "Job"))
        cmd = "\n".join(cmds)
        
        temdict = {
                "queue": self.pbs_queue,
                "walltime": self.pbs_walltime,
                "ppn": self.pbs_ppn,
                "email": self.pbs_email,
                "cmds": cmd
        }
        fo = "21.pbs"
        fho = open(fo, "w")
        assert op.isfile(self.pbs_template), "cannot read template: %s" % self.pbs_template
        fht = open(self.pbs_template, "r")
        src = Template(fht.read())
        fho.write(src.substitute(temdict))
        
        init()
        print(Fore.GREEN)
        print("A job script has been generated: %s" % fo)
        print("Please check, make necessary changes , then type:")
        print(Fore.RED + "qsub %s" % fo)
        print(Style.RESET_ALL)
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))
class shortread4Check(luigi.Task):
    name = luigi.Parameter(default="shortread4Check")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    ilist = luigi.Parameter(default="00.2.trim.tsv")
    olist = luigi.Parameter(default="00.3.tophat.tsv")
    def requires(self):
        return [shortread3Tophat(), WaitJob("%s/cps/shortread3TophatJob" % self.dirw)]
    def run(self):
        name, dirw, paired, ilist, olist = self.name, self.dirw, self.paired, self.ilist, self.olist
        os.chdir(dirw)
        assert op.isfile(ilist), "%s not exist" % ilist
        ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
        cols = list(ary.dtype.names)
        d22 = "22.tophat"
        d22 = op.abspath(d22)
        fho = open(olist, "w")
        if paired:
            fho.write("\t".join(cols + ["BamFile", "MappedPairs", "OrphanPairs", "UnmappedPairs", "UniquelyMappedPairs", "UniquelyMappedOrphans", "InsertSizeMean", "insertSizeStd"]) + "\n")
        else:
            fho.write("\t".join(cols + ["BamFile", "Mapped", "Unmapped", "UniquelyMapped"]) + "\n")
        for row in ary:
            row = [str(x, 'utf-8') for x in list(row)]
            sid = row[0]
            fbam = "%s/%s/accepted_hits.bam" % (d22, sid)
            assert check_bam(fbam), "%s not exist" % fbam
            if paired:
                f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
                fsum = "%s/%s/align_summary.txt" % (d22, sid)
            else:
                fr, rc, ft, rrc = row[5:9]
                fsum = "%s/%s/align_summary.txt" % (d22, sid)
                r1 = parse_tophat_alignsum(fsum)
                fsta = "%s/%s/samtools.stat" % (d22, sid)
                #r2 = parse_samtools_stat(fsta)
                assert rrc == r1['Input'], r1
                #assert r1['Mapped'] == r2['sequences'], (r1, r2)
                mapped = r1['Mapped']
                unmapped = str(int(rrc) - int(mapped))
                uni = str(int(mapped) - int(r1['of these']))
                fho.write("\t".join(row + [fbam, mapped, unmapped, uni])+"\n")
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))
class shortread5Htseq(luigi.Task):
    name = luigi.Parameter(default="shortread5Htseq")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    stranded = luigi.BoolParameter()
    ilist = luigi.Parameter(default="00.3.tophat.tsv")
    annotation = luigi.Parameter()
    pbs_template = luigi.Parameter()
    pbs_queue = luigi.Parameter()
    pbs_walltime = luigi.Parameter()
    pbs_ppn = luigi.Parameter()
    pbs_email = luigi.Parameter()
    def requires(self):
        return shortread4Check()
    def run(self):
        name, dirw, paired, stranded, ilist, annotation = self.name, self.dirw, self.paired, self.stranded, self.ilist, self.annotation
        os.chdir(dirw)
        assert op.isfile(ilist), "%s not exist" % ilist
        assert op.isfile(annotation), "%s not exist" % annotation
        ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
        fo1 = "31.1.htseq.sh"
        d32 = "32.htseq"
        fho1 = open(fo1, "w")
        for diro in [d32]:
            if not op.isdir(diro): 
                os.makedirs(diro)
        srd = 'no'
        if stranded: srd = 'yes'
        for row in ary:
            row = [str(x, 'utf-8') for x in list(row)]
            sid = row[0]
            if paired:
                fbam = row[15]
            else:
                fbam = row[9]
            fho1.write("samtools view %s | htseq-count -s %s -t gene -i ID -m union -a 20 - %s > %s/%s.txt\n" % (fbam, srd, annotation, d32, sid))
        
        cmds = []
        cmds.append("module load parallel")
        cmds.append("module load samtools")
        cmds.append("module load htseq/0.5.3") # requires python2
        cmds.append("cd %s" % dirw)
        cmds.append("parallel -j %s < %s" % (self.pbs_ppn, fo1))
        cmds.append("touch cps/%s" % (name + "Job"))
        cmd = "\n".join(cmds)
        
        temdict = {
                "queue": self.pbs_queue,
                "walltime": self.pbs_walltime,
                "ppn": self.pbs_ppn,
                "email": self.pbs_email,
                "cmds": cmd
        }
        fo = "31.pbs"
        fho = open(fo, "w")
        assert op.isfile(self.pbs_template), "cannot read template: %s" % self.pbs_template
        fht = open(self.pbs_template, "r")
        src = Template(fht.read())
        fho.write(src.substitute(temdict))
        
        init()
        print(Fore.GREEN)
        print("A job script has been generated: %s" % fo)
        print("Please check, make necessary changes , then type:")
        print(Fore.RED + "qsub %s" % fo)
        print(Style.RESET_ALL)
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))
class shortread6Check(luigi.Task):
    name = luigi.Parameter(default="shortread6Check")
    dirw = luigi.Parameter()
    paired = luigi.BoolParameter()
    ilist = luigi.Parameter(default="00.3.tophat.tsv")
    olist = luigi.Parameter(default="00.5.htseq.tsv")
    def requires(self):
        return [shortread5Htseq(), WaitJob("%s/cps/shortread5HtseqJob" % self.dirw)]
    def run(self):
        name, dirw, paired, ilist, olist = self.name, self.dirw, self.paired, self.ilist, self.olist
        os.chdir(dirw)
        assert op.isfile(ilist), "%s not exist" % ilist
        ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
        cols = list(ary.dtype.names)
        d32 = op.abspath("32.htseq")
        fho = open(olist, "w")
        fho.write("\t".join(cols + ["HtseqFile"])+"\n")
        for row in ary:
            row = [str(x, 'utf-8') for x in list(row)]
            sid = row[0]
            fhts = "%s/%s.txt" % (d32, sid)
            assert op.isfile(fhts), "%s not there" % fhts
            if paired:
                f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            else:
                fr, rc, ft, rrc = row[5:9]
            fho.write("\t".join(row + [fhts]) + "\n")
        os.system("touch %s/cps/%s" % (dirw, name))
    def output(self):
        return luigi.LocalTarget("%s/cps/%s" % (self.dirw, self.name))

if __name__ == "__main__":
    luigi.run()

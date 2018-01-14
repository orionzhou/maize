#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import time
import argparse

def run_cmds_ssh(cmds):
    username = 'zhoux379'
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
            print("... " + line.strip("\n"))
    s.close()

def genome_fas(org, rename):
    dirw = op.join("/home/springer/zhoux379/data/genome", org)
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    for fname in ["raw.fix.fas.index", "11_genome.fas.index"]:
        if op.isfile(fname):
            os.remove(fname)
    if op.islink("11_genome.fas"): os.unlink("11_genome.fas")
        
    assert op.isfile("raw.fas"), "no raw.fas in %s" % dirw
    os.system("seq.check.pl -i raw.fas -o raw.fix.fas")
    if rename:
        os.system("seq.rename.pl -i raw.fix.fas -p scf -o 11_genome.fas")
        os.system("rm raw.fix.fas")
    else:
        os.system("mv raw.fix.fas 11_genome.fas")
    if op.isfile("ctg.raw.fas"):
        os.system("seq.check.pl -i ctg.raw.fas -o ctg.fas")
    os.system("seqlen.py 11_genome.fas 15.sizes")
    os.system("awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print $1, 0, $2}' 15.sizes > 15.bed")
    os.system("seqgap.pl -i 11_genome.fas -o 16.gap.bed -m 10")
    #os.system("bedToBigBed 16.gap.bed 15.sizes 16.gap.bb")
    #os.system("bgzip -c 16.gap.bed > 16.gap.bed.gz")
    #os.system("tabix -p bed 16.gap.bed.gz")
        
    if op.isdir("12.rm"): os.system("rm -rf 12.rm")
    os.makedirs("12.rm")
    
def build_blat(org):
    dirg = op.join("/home/springer/zhoux379/data/genome", org)
    fg = "%s/11_genome.fas" % dirg
    dirw = op.join(dirg, "21.blat")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    
    os.system("faToTwoBit %s db.2bit" % fg)
    os.system("blat db.2bit tmp.fas tmp.out -makeOoc=db.2bit.tile11.ooc")
    if op.isfile("tmp.out"): os.remove("tmp.out")
def build_bowtie2(org):
    dirg = op.join("/home/springer/zhoux379/data/genome", org)
    fg = "%s/11_genome.fas" % dirg
    dirw = op.join(dirg, "21.bowtie2")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    
    if op.isfile("db.fas"): os.system("rm db.fas")
    os.system("ln -sf %s db.fas" % fg)
    # need to "module load bowtie2"
    os.system("bowtie2-build db.fas db")
def build_bwa(org):
    dirg = op.join("/home/springer/zhoux379/data/genome", org)
    fg = "%s/11_genome.fas" % dirg
    dirw = op.join(dirg, "21.bwa")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    os.system("bwa index -p %s/db %s" % (dirw, fg))
def repeatmasker_parse(org):
    dirw = op.join(os.environ['genome'], org)
    os.chdir(dirw)
    os.system("parse.rm.pl -i 12.rm/11_genome.fas.out -o 12.rm.tbl")
    os.system("awk 'BEGIN{OFS=\"\\t\"} {print $1, $2-1, $3, $9 \" | \" ($5)}' 12.rm.tbl > 12.rm.raw.bed")
    os.system("sortBed -i 12.rm.raw.bed > 12.rm.bed")
    #os.system("bedToBigBed -tab -type=bed4 12.rm.bed 15.sizes 12.rm.bb")
    os.system("rm 12.rm.raw.bed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'process genome files and build genome DB'
    )
    parser.add_argument(
            'org', help = 'species name / prefix'
    )
    parser.add_argument(
            '--fasta', action='store_true', help = 'process Fasta [No]'
    )
    parser.add_argument(
            '--rename', action='store_true', help = 'rename scaffold IDs [No]'
    )
    parser.add_argument(
            '--blat', action='store_true', help='built BLAT DB [No]'
    )
    parser.add_argument(
            '--bowtie2', action='store_true', help='built Bowtie2 DB [No]'
    )
    parser.add_argument(
            '--bwa', action='store_true', help='built BWA DB [No]'
    )
    args = parser.parse_args()
    org = args.org

    dirg = op.join("/home/springer/zhoux379/data/genome", org)
    dird = "/home/springer/zhoux379/data/db"
    if args.fasta: genome_fas(org, args.rename)
    if args.blat: build_blat(org)
    if args.bowtie2: build_bowtie2(org)
    if args.bwa: build_bwa(org)

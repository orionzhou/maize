#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import time
import logging

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

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

def clean_fasta(args):
    dirw = op.join("/home/springer/zhoux379/data/genome", args.species)
    if not op.isdir(dirw): makedir(dirw)
    os.chdir(dirw)
    for fname in ["raw.fix.fas.index", "11_genome.fas.index"]:
        if op.isfile(fname):
            os.remove(fname)
    if op.islink("11_genome.fas"): os.unlink("11_genome.fas")
   
    if op.isfile("11_genome.fas") and not args.overwrite:
        logging.debug("11_genome.fas already exits: skipped")
    else:
        fis = ['raw.fas', 'raw.fa', 'raw.fas.gz', 'raw.fa.gz']
        fis = [x for x in fis if op.isfile(x)]
        if len(fis) == 0:
            eprint("no raw.fas found")
            sys.exit(1)
        elif len(fis) > 1:
            eprint(">1 raw.fas found")
            sys.exit(1)
        fi = fis[0]
        if fi.endswith(".gz"):
            sh("gunzip -c %s | fasta clean - > 01.fas" % fi)
        else:
            sh("fasta clean %s > 01.fas" % fi)

        if not args.norename:
            sh("fasta rename --map 03.seqid.map 01.fas > 11_genome.fas")
            os.remove("01.fas")
        else:
            sh("mv 01.fas 11_genome.fas")
    
    if op.isfile("ctg.raw.fas"):
        sh("fasta clean ctg.raw.fas > ctg.fas")

    if op.isfile("15.sizes") and not args.overwrite:
        logging.debug("15.sizes already exits - skipped")
    else:
        sh("fasta size 11_genome.fas > 15.sizes")
    
    if op.isfile("15.bed") and not args.overwrite:
        logging.debug("15.bed already exits - skipped")
    else:
        sh("fasta size --bed 11_genome.fas > 15.bed")
    
    if op.isfile("16.gap.bed") and not args.overwrite:
        logging.debug("16.gap.bed already exits - skipped")
    else:
        sh("fasta gaps 11_genome.fas > 16.gap.bed")

def build_blat(args):
    dirg = op.join("/home/springer/zhoux379/data/genome", args.species)
    fg = "%s/11_genome.fas" % dirg
    assert op.isfile(fg), "%s not there" % fg
    dirw = op.join(dirg, "21.blat")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
   
    if not args.overwrite and op.isfile('db.2bit'):
        logging.debug("db.2bit already exists - skipped")
    else:
        sh("faToTwoBit %s db.2bit" % fg)
        sh("blat db.2bit tmp.fas tmp.out -makeOoc=db.2bit.tile11.ooc")
    if op.isfile("tmp.out"): os.remove("tmp.out")

def build_bowtie(args):
    dirg = op.join("/home/springer/zhoux379/data/genome", args.species)
    fg = "%s/11_genome.fas" % dirg
    assert op.isfile(fg), "%s not there" % fg
    dirw = op.join(dirg, "21.bowtie2")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    
    if op.isfile("db.rev.1.bt2") and not args.overwrite:
        logging.debug("db.*.bt2 already exists - skipped")
    else:
        sh("rm -rf *")
        sh("ln -sf %s db.fas" % fg)
        # need to "module load bowtie2"
        sh("bowtie2-build db.fas db")

def build_hisat(args):
    dirg = op.join("/home/springer/zhoux379/data/genome", args.species)
    fg = "%s/11_genome.fas" % dirg
    assert op.isfile(fg), "%s not there" % fg
    dirw = op.join(dirg, "21.hisat")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
   
    if op.isfile("db.bwt") and not args.overwrite:
        logging.debug("db.bwt already exists - skipped")
    else:
        sh("bwa index -p %s/db %s" % (dirw, fg))

def build_bwa(args):
    dirg = op.join("/home/springer/zhoux379/data/genome", args.species)
    fg = "%s/11_genome.fas" % dirg
    assert op.isfile(fg), "%s not there" % fg
    dirw = op.join(dirg, "21.bwa")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
   
    if op.isfile("db.bwt") and not args.overwrite:
        logging.debug("db.bwt already exists - skipped")
    else:
        sh("bwa index -p %s/db %s" % (dirw, fg))

def repeatmasker(args):
    from maize.formats.pbs import PbsJob
    dirg = op.join("/home/springer/zhoux379/data/genome", args.species)
    fg = "%s/11_genome.fas" % dirg
    dirw = op.join(dirg, "12.repeatmasker")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)

    species = None
    if args.species in ['Zmays', 'B73', 'PH207', 'W22', 'Mo17']:
        species = 'maize'
    elif args.species == 'Osativa':
        species = 'rice'
    else:
        logging.error("%s not supported" % args.species)
        sys.exit(1)
    cmds = [
        "cd %s" % dirw,
        "RepeatMasker -pa %d -species %s -dir %s %s" % (args.cpu, species, dirw, fg),
        "parse.rm.pl -i 11_genome.fas.out -o 12.repeatmasker.tsv"
    ]
    pbsjob = PbsJob(ppn = 24, walltime = "5:00:00", cmds = "\n".join(cmds))
    fjob = op.join(dirg, "13.rm.pbs")
    pbsjob.write(fjob)
    logging.debug("Job script '%s' has been created" % fjob)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'process genome files and build genome DB'
    )
    parser.add_argument('species', help = 'species/accession/genotype name')
    parser.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("fasta", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "clean and rename fasta records, generate *.sizes and gap location files")
    sp1.add_argument('--norename', action = 'store_true', help = 'don\'t rename seq IDs')
    sp1.set_defaults(func = clean_fasta)

    sp1 = sp.add_parser("repeatmasker", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "run repeatmasker and parse result"
    )
    sp1.add_argument('--cpu', default = 24, help = 'number CPUs to use')
    sp1.set_defaults(func = repeatmasker)

    sp2 = sp.add_parser("blat", help = "build Blat DB")
    sp2.set_defaults(func = build_blat)
    sp2 = sp.add_parser("bowtie", help = "build Bowtie2 DB")
    sp2.set_defaults(func = build_bowtie)
    sp2 = sp.add_parser("bwa", help = "build bwa DB")
    sp2.set_defaults(func = build_bwa)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()
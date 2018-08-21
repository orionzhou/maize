#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import time
import logging

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob

def make_genomedir(species):
    dirw = species
    if species.isalnum():
        dirw = op.join("/home/springer/zhoux379/data/genome", species)
        logging.debug("converting species to directory: %s" % dirw)
    if not op.isdir(dirw):
        logging.debug("creating diretory: %s" % dirw)
        mkdir(dirw)
    return dirw

def get_genomedir(species):
    dirw = species
    if species.isalnum():
        dirw = op.join("/home/springer/zhoux379/data/genome", species)
        logging.debug("converting species to directory: %s" % dirw)
    fg = "%s/10_genome.fna" % dirw
    if not op.isfile(fg):
        logging.error("%s not there" % fg)
        sys.exit()
    return op.abspath(dirw), op.abspath(fg)

def clean_fasta(args):
    dirw = make_genomedir(args.species)
    os.chdir(dirw)
    for fname in ["raw.fix.fas.index", "11_genome.fas.index"]:
        if op.isfile(fname):
            os.remove(fname)
    if op.islink("10_genome.fna"): os.unlink("10_genome.fna")
   
    if op.isfile("10_genome.fna") and not args.overwrite:
        logging.debug("10_genome.fna already exits: skipped")
    elif op.isfile("08_seq_map/renamed.fna"):
        sh("ln -sf 08_seq_map/renamed.fna 10_genome.fna")
        if op.isfile("08_seq_map/renamed.sizes"):
            sh("ln -sf 08_seq_map/renamed.sizes 10_genome.sizes")
    else:
        logging.error("08_seq_map/renamed.fna not there")
        sys.exit(1)
    
    if not op.isdir("15_intervals"):
        mkdir("15_intervals")
    
    if op.isfile("15_intervals/01.chrom.bed") and not args.overwrite:
        logging.debug("01.chrom.bed already exits - skipped")
    else:
        sh("fasta size --bed 10_genome.fna > 15_intervals/01.chrom.bed")
    if op.isfile("15_intervals/01.chrom.sizes") and not args.overwrite:
        logging.debug("01.chrom.sizes already exits - skipped")
    else:
        sh("faSize -detailed 10_genome.fna > 15_intervals/01.chrom.sizes")
    
    if op.isfile("15_intervals/11.gap.bed") and not args.overwrite:
        logging.debug("11.gap.bed already exits - skipped")
    else:
        sh("fasta gaps 10_genome.fna > 15_intervals/11.gap.bed")

def build_blat(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "21_dbs/blat")
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
   
    if not args.overwrite and op.isfile('db.2bit'):
        logging.debug("db.2bit already exists - skipped")
    else:
        sh("faToTwoBit %s db.2bit" % fg)
        sh("blat db.2bit tmp.fas tmp.out -makeOoc=db.2bit.tile11.ooc")
    if op.isfile("tmp.out"): os.remove("tmp.out")

def build_bowtie(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "21_dbs/bowtie2")
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
    
    if op.isfile("db.rev.1.bt2") and not args.overwrite:
        logging.debug("db.*.bt2 already exists - skipped")
    else:
        sh("rm -rf *")
        sh("ln -sf %s db.fa" % fg)
        # need to "module load bowtie2"
        sh("bowtie2-build db.fa db")

def build_hisat(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "21_dbs/hisat2")
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
   
    f_gtf = "../../50_annotation/10.gtf"
    if op.isfile("db.1.ht2") and not args.overwrite:
        logging.debug("db.1.ht2 already exists - skipped")
    elif not op.isfile(f_gtf):
        logging.error("no gtf file: f_gtf")
        sys.exit()
    else:
        sh("hisat2_extract_exons.py %s > db.exon" % f_gtf)
        sh("hisat2_extract_splice_sites.py %s > db.ss" % f_gtf)
        sh("hisat2-build -p %d --ss db.ss --exon db.exon %s db" % (args.p, fg))

def build_star(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "21_dbs/star")
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
  
    f_gtf = "../../50_annotation/10.gtf"
    if op.isfile("SA") and not args.overwrite:
        logging.debug("SA already exists - skipped")
    elif not op.isfile(f_gtf):
        logging.error("no gtf file: %s" % f_gtf )
        sys.exit()
    else:
        sh("STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s \
                --genomeFastaFiles %s --sjdbGTFfile %s" %
                (args.p, ".", fg, f_gtf))

def build_bwa(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "21_dbs/bwa")
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
   
    if op.isfile("db.bwt") and not args.overwrite:
        logging.debug("db.bwt already exists - skipped")
    else:
        sh("bwa index -a bwtsw -p %s/db %s" % (dirw, fg))

def build_gatk(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "21_dbs/gatk")
    if not op.isdir(dirw): mkdir(dirw)
    os.chdir(dirw)
   
    if op.isfile("db.dict") and not args.overwrite:
        logging.debug("db.dict already exists - skipped")
    else:
        if op.exists("db.fasta"): sh("rm db.fasta")
        if op.exists("db.dict"): sh("rm db.dict")
        sh("cp ../../10_genome.fna db.fasta")
        sh("gatk CreateSequenceDictionary -R db.fasta")
        sh("samtools faidx db.fasta")
        #sh("gatk FindBadGenomicKmersSpark -R db.fasta -O kmers_to_ignore.txt")
        #sh("gatk BwaMemIndexImageCreator -I db.fasta -O db.img")

def repeatmasker(args):
    dirg, fg = get_genomedir(args.species)
    dirw = op.join(dirg, "12_repeatmasker")
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)

    species = None
    if args.species in ['Zmays', 'B73', 'PH207', 'W22', 'Mo17', 'PHB47']:
        species = 'maize'
    elif args.species == 'Osativa':
        species = 'rice'
    else:
        logging.error("%s not supported" % args.species)
        sys.exit(1)
    
    cmds = []
    cmds.append("cd %s" % dirw)
    cmds.append("RepeatMasker -pa %d -species %s -dir %s %s" % (args.p, species, dirw, fg)),
    cmds.append("parse.rm.pl -i 11_genome.fas.out -o 12.repeatmasker.tsv")
    
    pbsjob = PbsJob(queue = 'ram256g', ppn = 24, walltime = "10:00:00", cmds = "\n".join(cmds))
    fjob = op.join(dirg, "13.rm.pbs")
    pbsjob.write(fjob)
    logging.debug("Job script '%s' has been created" % fjob)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'process genome files and build genome DB'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("fasta", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "clean and rename fasta records, generate *.sizes and gap location files")
    sp1.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp1.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp1.set_defaults(func = clean_fasta)

    sp1 = sp.add_parser("repeatmasker", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "run repeatmasker and parse result"
    )
    sp1.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp1.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp1.add_argument('--p', type = int, default = 24, help = 'number of threads')
    sp1.set_defaults(func = repeatmasker)

    sp2 = sp.add_parser("blat", help = "build Blat DB")
    sp2.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp2.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp2.set_defaults(func = build_blat)
    
    sp2 = sp.add_parser("bowtie", help = "build Bowtie2 DB")
    sp2.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp2.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp2.set_defaults(func = build_bowtie)
    
    sp2 = sp.add_parser("bwa", help = "build bwa DB")
    sp2.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp2.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp2.set_defaults(func = build_bwa)
    
    sp2 = sp.add_parser("hisat", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "build hisat2 DB"
    )
    sp2.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp2.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp2.add_argument('--p', type = int, default = 1, help = 'number of threads')
    sp2.set_defaults(func = build_hisat)
    
    sp2 = sp.add_parser("star", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "build STAR DB"
    )
    sp2.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp2.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp2.add_argument('--p', type = int, default = 1, help = 'number of threads')
    sp2.set_defaults(func = build_star)
    
    sp2 = sp.add_parser("gatk", help = "build GATK ref-db")
    sp2.add_argument('species', help = 'species/accession/genotype/dir-path')
    sp2.add_argument('--overwrite', action='store_true', help = 'overwrite')
    sp2.set_defaults(func = build_gatk)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

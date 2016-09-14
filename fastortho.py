#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_cfg(fc):
    (orgs, fis) = ([], [])
    fhc = open(fc, "r")
    for line in fhc:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        (org, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            print os.access(fi, os.F_OK)
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
    fhc.close()
    return (orgs, fis)
def merge_seqs(fc, fo):
    print "merging input files to %s" % fo
    (orgs, fis) = read_cfg(fc)
    seqs = []
    for i in range(0,len(orgs)):
        handle = 0
        if (fis[i].endswith(".gz")):
            handle = gzip.open(fis[i], "rb")
        else:
            handle = open(fis[i], "rU")
        seq_it = SeqIO.parse(handle, "fasta")
        handle.close

        seqs1 = [SeqRecord(rcd.seq, id = 
            ("%s|%s [%s.]" % (orgs[i], rcd.id, orgs[i])),
            description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()
def write_config(fo, f_fas, f_blast):
    fho = open(fo, "w")
    print >>fho, "--mixed_genome_fasta %s" % f_fas
    print >>fho, "--working_directory ."
    print >>fho, "--project_name zz"
    print >>fho, "--blast_file %s" % f_blast
    fho.close()
def reformat_blast_id(fi, fo):
    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        line = line.strip("\n")
        line = line.strip("\r")
        ary = line.split("\t")
        (org1, gid1) = ary[0].split("|")
        (org2, gid2) = ary[1].split("|")
        ary[0] = "%s [%s]" % (gid1, org1)
        ary[1] = "%s [%s]" % (gid2, org2)
        print >>fho, "\t".join(ary) 
    fhi.close()
    fho.close()
def fastortho2tbl(fi, fs, fo):
    fhs = open(fs, "r")
    aids = set()
    for line in fhs:
        line = line.strip("\n")
        ary = line.split("\t")
        aids.add(ary[0])
    fhs.close()

    fhi = open(fi, "r")
    fho = open(fo, "w")
    i = 1
    print >>fho, "grp\tid"
    for line in fhi:
        line = line.strip("\n")
        line = line.strip("\r")
        ary = line.split(":")
        ary = ary[1].strip().split(" ")
        for ele in ary:
            (id, tmp) = ele.split("(")
            print >>fho, "%d\t%s" % (i, id)
            aids.remove(id)
        i += 1
    fhi.close()

    for rid in aids:
        print >>fho, "%d\t%s" % (i, rid)
        i += 1
    fho.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        'cfgfile', help = 'config file (.csv)'
    )
    parser.add_argument(
        'outdir', help = 'output directory'
    )
    args = parser.parse_args()
    (fc, dirw) = (args.cfgfile, args.outdir)
    dirw = op.realpath(dirw)
    
    if not op.exists(dirw): os.makedirs(dirw)
    
    cdir = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] = os.environ['PATH']+':'+cdir
    cwd = os.getcwd()
    os.chdir(dirw)
    
    #merge_seqs(fc, "01.fas")
    #os.system("seqlen.py 01.fas 02.tbl")
    #os.system("blast.self.py 01.fas 06 --cpu 8")
    
    #write_config("config", "01.fas", '06.out')
    os.system("FastOrtho --option_file config")
    fastortho2tbl("zz.end", "02.tbl", "11.tbl")
   

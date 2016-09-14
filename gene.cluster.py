#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
import os.path as op
import numpy as np
import argparse

orgs = """HM101
    HM058 HM125 HM056 HM129 HM060
    HM095 HM185 HM034 HM004 HM050
    HM023 HM010 HM022 HM324 HM340
    HM056.AC HM034.AC HM340.AC""".split()
dirw = "/home/youngn/zhoup/Data/misc2/gene.cluster"
if not op.exists(dirw): os.makedirs(dirw)

def make_blastp_db(orgs):
    dirg = os.environ['genome']
    diro = op.join(os.environ['data'], "db", "blastp")
    if not op.exists(diro): os.makedirs(diro)
   
    for org in orgs:
        ff = op.join(dirg, org, "51.fas")
        fo = op.join(diro, org)
        cmd = "makeblastdb -dbtype prot -in %s -out %s" % (ff, fo)
        #print cmd
        os.system(cmd)
    return
def run_blast(org, fo):
    dirg = op.join(os.environ['genome'], org)
    ff = op.join(dirg, "51.fas")
    fdb = op.join(os.environ['data'], 'db', 'blastp', org)
    cmd = "blastp -db %s -query %s -out %s -num_threads %d -evalue 1e-5 -outfmt '7 qseqid qlen sseqid slen length bitscore evalue'" % (fdb, ff, fo, 24)
    print cmd
def blast2tbl(fi, fo):
    e_min = 1
    fhi = open(fi, "r")
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, qlen, tid, tlen, alen, bit, e) = line.split("\t")
        e = float(e)
        if e != 0 and e < e_min: e_min = e
    fhi.close()
    
    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, qlen, tid, tlen, alen, bit, e) = line.split("\t")
        (qlen, tlen, alen) = (float(qlen), float(tlen), float(alen))
        (e, bit) = (float(e), float(bit))
        if e == 0: e = e_min
        rlen = alen / min(qlen, tlen)
        if qid != tid and rlen >= 0.5 and e < 1e-10:
            score = - math.log10(e)
            print >>fho, "%s\t%s\t%g" % (qid, tid, score)
    fhi.close()
    fho.close()
def run_blat(dirw, org):
    dirg = os.environ['genome']
    diro = op.join(dirw, "01_pro_blat")
    if not op.exists(diro): os.makedirs(diro)
    
    ff = op.join(dirg, org, "51.fas")
    fo = op.join(diro, org + ".tbl")
    cmd = "blat -prot -out=blast8 %s %s %s" % (ff, ff, fo)
    print cmd
    os.system(cmd)
    return
def blat_filter(dirw, org):
    diri = op.join(dirw, "01_pro_blat")
    diro = op.join(dirw, "02_filtered")
    if not op.exists(diro): os.makedirs(diro)
    
    fi = op.join(diri, org + ".tbl")
    fo = op.join(diro, org + ".tbl")
    
    fhi = open(fi, "r")
    fho = open(fo, "w")
    
    for line in fhi:
        line = line.strip("\n")
        (qid, tid, idt, alen, mis, gap, qbeg, qend, tbeg, tend, e, bit) = line.split("\t")
        (e, bit) = (float(e), float(bit))
        if qid != tid and e < 1e-5 and alen:
            print >>fho, line
    fhi.close()
    fho.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='JJJ')
    parser.add_argument('--version', action="version", version="%(prog)s 1.0")
    args = parser.parse_args()
    
    d01 = op.join(dirw, "01.blast")
    d03 = op.join(dirw, "03.mcl.in")
    d04 = op.join(dirw, "04.mcl")
    d08 = op.join(dirw, "08.grp")
    for dirw in [d01, d03, d04, d08]:
        if not op.exists(dirw): os.makedirs(dirw)
    
    for org in orgs:
        print "working on " + org
        f01 = op.join(d01, org + ".tbl")
        #run_blast(org, f01)

        f03 = op.join(d03, org + ".tbl")
        blast2tbl(f01, f03)

        f04 = op.join(d04, org + ".mcl")
        cmd = "$soft/mcl/bin/mcl %s -te 4 -I 2.0 --abc -o %s" % (f03, f04)
        os.system(cmd)
        
        f08 = op.join(d08, org + ".tbl")

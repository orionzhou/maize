#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
import re 
import os.path as op
import numpy as np
import argparse

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
    (cqid, ctid, cqlen, ctlen, calen, cscore) = ('', '', 0, 0, 0, 0)
    for line in fhi:
        if line[0] == "#":
            continue
        line = line.strip("\n")
        (qid, qlen, tid, tlen, alen, bit, e) = line.split("\t")
        if qid == tid: continue
        (qlen, tlen, alen) = (float(qlen), float(tlen), float(alen))
        (e, bit) = (float(e), float(bit))
        if e == 0: e = e_min
        score = - math.log10(e)
        if qid == cqid and tid == ctid:
            calen += alen
            cscore += score / 2
        else:
            if cqid != '':
                if calen/min(cqlen, ctlen) >= 0.5 and cscore > 10:
                    print >>fho, "%s\t%s\t%g" % (cqid, ctid, cscore)
            (cqid, ctid, cqlen, ctlen, calen, cscore) = (qid, tid, qlen, tlen, alen, score)
    if calen/min(cqlen, ctlen) >= 0.5 and cscore > 10:
        print >>fho, "%s\t%s\t%g" % (cqid, ctid, cscore)
    fhi.close()
    fho.close()
def mcl2tbl(fi, fo):
    fhi = open(fi, "r")
    fho = open(fo, "w")
    print >>fho, "grp\tid"
    clu = 1
    for line in fhi:
        line = line.strip("\n")
        lst = line.split("\t")
        for ele in lst:
            print >>fho, "%d\t%s" % (clu, ele)
        clu += 1
    fhi.close()
    fho.close()
def cdhit2tbl(fi, fo):
    fhi = open(fi, "r")
    fho = open(fo, "w")
    print >>fho, "grp\tid"

    grp = 0
    p = re.compile('\d+.*>(\S+)\.{3}\s')
    for line in fhi:
        if line[0] == ">":
            grp += 1
            continue
        line = line.strip("\n")
        m = p.match(line)
        print >>fho, "%d\t%s" % (grp, m.group(1))
    fhi.close()
    fho.close()
def usearch2tbl(fi, fo):
    fhi = open(fi, "r")
    fho = open(fo, "w")
    fho.write("grp\tid\n")

    for line in fhi:
        line = line.strip("\n")
        (tag, grp, len_size, idty, srd, ign1, ign2, cigar, id, tgt) = line.split("\t")
        grp = int(grp) + 1
        if tag == "H" or tag == "S":
            fho.write("%d\t%s\n" % (grp, id))
    fhi.close()
    fho.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'cluster a set of protein sequences'
    )
    parser.add_argument(
        'fi', help = 'input file (fasta)'
    )
    parser.add_argument(
        'dirw', help = 'output directory'
    )
    parser.add_argument(
        'idty', nargs = '?', help = 'clustering similarity (0.7)'
    )
    args = parser.parse_args()
   
    idty = 0.7
    (fi, dirw) = (args.fi, args.dirw)
    if args.idty: idty = float(args.idty)
    if not op.exists(dirw): os.makedirs(dirw)
    
    f01 = op.join(dirw, "01.blastdb")
    cmd = "makeblastdb -dbtype prot -in %s -out %s" % (fi, f01)
#    os.system(cmd)
   
    f11 = op.join(dirw, "11.blast.tbl")
    cmd = "blastp -db %s -query %s -out %s -num_threads %d -evalue 1e-5 -outfmt '7 qseqid qlen sseqid slen length bitscore evalue'" % (f01, fi, f11, 24)
#    print cmd
#    os.system(cmd)
    
    f12 = op.join(dirw, "12.mcl.in.tbl")
#    blast2tbl(f11, f12)

    f15 = op.join(dirw, "15.mcl")
    cmd = "$soft/mcl/bin/mcl %s -te 4 -I 8.0 --abc -o %s" % (f12, f15)
#    os.system(cmd)
    
    f16 = op.join(dirw, "16.tbl")
#    mcl2tbl(f15, f16)
    
    f21 = op.join(dirw, "21.cdhit")
    f22 = op.join(dirw, "22.tbl")
    cmd = "$git/cdhit/cd-hit -c %g -T 0 -aL 0.5 -n 2 -M 0 -d 0 -i %s -o %s" % (idty, fi, f21)
#    os.system(cmd)
#    cdhit2tbl(f21 + ".clstr", f22)

    f31 = op.join(dirw, "31.uc")
    f32 = op.join(dirw, "32.tbl")
    cmd = "usearch -cluster_fast %s -sort length -id %g -uc %s" % (fi, idty, f31)
    os.system(cmd)
    usearch2tbl(f31, f32)

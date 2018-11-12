#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord

def parse_gi(fi):
    fhi = open(fi, "r")
    res = dict()
    for line in fhi:
        line = line.strip("\n")
        gi = line.split("|")[1]
        if res.has_key(gi):
            print("  %s already appeared - skipped" % gi)
            continue
        res[gi] = []
        tmp = line.split(":")
        if len(tmp) > 1:
            beg, end = map(int, tmp[1].split("-"))
            res[gi] = [beg, end]
    fhi.close()
    print("%d unique GIs extracted" % len(res.keys()))
    return res

def gi_seqret(res, fs, fo):
    fho = open(fo, "w")

    Entrez.email="zhoupenggeni@gmail.com"
    seqs = []
    for gi in res.keys():
        loc = res[gi]
        handle = Entrez.efetch(db="nucleotide", id=gi, rettype="fasta")
        orcd = SeqIO.read(handle, "fasta")
        handle.close()
        if len(loc) == 2:
            beg, end = loc
        else:
            beg, end = [1, len(orcd)]
        sid = "%s-%d-%d" % (gi, beg, end)
        rcd = SeqRecord(orcd.seq[beg-1:end], id=sid, description='')
        seqs.append(rcd)

        desc = orcd.description
        eles = desc.split(" ")
        name1, name2 = eles[1:3]
        print >>fho, "%s\t%d\t%d\t%s\t%s\t%s" % (gi, beg, end, name1, name2, desc)
    fho.close()

    fhs = open(fs, "w")
    SeqIO.write(seqs, fhs, "fasta")
    fhs.close()

dirw = '/home/youngn/zhoup/Data/misc2/bact'
if __name__ == "__main__":
    res = parse_gi("gis.txt")
    gi_seqret(res, "11.fas", "11.tsv")


#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
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
def internet_on():
    try:
        resp = urllib2.urlopen("http://www.google.com/", timeout = 3)
        return True
    except urllib2.URLError as err: pass
    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'blast a set of proteins against NR'
    )
    parser.add_argument(
        'fi', help = 'input file (fasta)'
    )
    parser.add_argument(
        'fo', help = 'output file (prefix)'
    )
    args = parser.parse_args()
    (fi, fo) = (args.fi, args.fo)
   
    nr = os.environ['data'] + "/db/blast/current/nr";
    nproc = int(os.environ['nproc'])

    f01 = fo + ".1.blast.tbl"
    cmd = """blastp -db %s -outfmt \
        '6 qseqid qstart qend qlen sseqid sstart send slen length nident mismatch gaps evalue bitscore' \
        -max_target_seqs 3 -evalue 1e-5 -num_threads %d \
        -query %s -out %s""" % (nr, nproc, fi, f01)
    print cmd
    os.system(cmd)

    if internet_on:
        print "internet on :-)"
    else:
        print "no internet :-("

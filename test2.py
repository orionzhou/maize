#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path as op
#from time import clock, time
from subprocess import Popen, PIPE
from Bio import SeqIO

def seqret(f_seq, id, beg, end, strand=1):
    if not os.path.isfile(f_seq):
        sys.exit(f_seq + " is not there")
    f_idx = f_seq + ".idx"
    seq_idx = SeqIO.index_db(f_idx, f_seq, "fasta")
    if not id in seq_idx:
        sys.exit("cannot find " + id + "in " + f_seq)
    
    seq_obj = seq_idx[id].seq[beg-1:end]
    if strand == -1 or strand == "-":
        seq_obj = seq_obj.reverse_complement()
    seq_str = seq_obj.tostring()
    return seq_str
def main2():
    seq_db = "/home/youngn/zhoux379/data/misc3/spada.crp/Athaliana/01_genome/01_refseq.fa"
    seq_id = "Chr3"
    itv = 100000
    seq_beg = itv * 1 + 1
    seq_end = itv * 1 + 10
        
    c0, t0 = clock(), time()
    print(seqret(seq_db, seq_id, seq_beg, seq_end, 1))

    print("process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0))

    c0, t0 = clock(), time()
    f_perl = "/soft/bioperl/5.16.1/bin/bp_seqret.pl"
    cmd = '%s -d %s -i %s -s %d -e %d' % (f_perl, seq_db, seq_id, seq_beg, seq_end)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    lines = p.stdout.readlines()
    print(lines)
    print("process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0))
def corncyc2tsv():
    dirw = '/home/springer/zhoux379/data/genome/Zmays_v4/corncyc'
    fi = op.join(dirw, "query-results.tsv")
    fo = op.join(dirw, "01.tsv")
    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        ps = line.strip().split("\t")
        if len(ps) < 2:
            continue
        else:
            name, gidstr = ps[0:2]
            name = name.strip("\"").replace("&","").replace(";","")
            name = name.replace("<i>","").replace("</i>","")
            name = name.replace("<sub>","").replace("</sub>","")
            gids = gidstr.strip("()").split(" ")
            gids = [x.strip("\"").split("_")[0] for x in gids]
            if len(gids) > 0:
                for gid in gids:
                    fho.write("%s\t%s\n" % (name, gid))
    fhi.close()
    fho.close()

if __name__ == "__main__":
    corncyc2tsv()

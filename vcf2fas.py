#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'make an alternate fasta file using vcf'
    )
    parser.add_argument(
        'fi', help = 'interval file (.tsv)'
    )
    parser.add_argument(
        'fv', help = 'vcf file (.vcf)'
    )
    parser.add_argument(
        'fs', help = 'reference sequence file (.fas)'
    )
    parser.add_argument(
        'fo', help = 'output file (.tsv)'
    )
    args = parser.parse_args()
    fi, fv, fs, fo = args.fi, args.fv, args.fs, args.fo
    
    fhi = open(fi, "r")
    fho = open(fo, "w")
    seqdic = SeqIO.index(fs, "fasta")
    vcfr = vcf.Reader(filename = fv, compressed = True)
     
    for line in fhi:
        line = line.strip("\n")
        ps = line.split("\t")
        if len(ps) < 4:
            print("<4 fields: %s" % line)
            continue
        seqid, beg, end, note = ps[0:4]
        beg, end = int(beg), int(end)
        chrstr = seqdic[seqid].seq
        seqstr = chrstr[beg-1:end]
        rcds = vcfr.fetch(seqid, start = beg-1, end = end)
        fho.write("%s\t%s\t%s\t%s\t%s\n" % (seqid, beg, end, note, seqstr))
        print("%s\t%s\t%s\t%s" % (seqid, beg, end, note))
        for rcd in rcds:
            alts = ",".join(map(str, rcd.ALT))
            alt = rcd.ALT[0]
            print("%s\t%s\t%s\t%s" % (rcd.CHROM, rcd.POS, rcd.REF, alts))

#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'translate mo17 cds sequences'
    )
    parser.add_argument(
            'fi', help = 'input file (.tsv)'
    )
    parser.add_argument(
            'fo', help = 'output file (.tsv)'
    )
    parser.add_argument(
            'fs', help = 'output fasta file (.fasta)'
    )
    args = parser.parse_args()
    fi, fo, fs = args.fi, args.fo, args.fs
    fhi, fho, fhs = open(fi, "r"), open(fo, "w"), open(fs, "w")
    
    rcds = []
    for line in fhi:
        line = line.strip("\n")
        ary = line.split("\t")
        if len(ary) != 4:
            print("not 4 fields:\n$line")
            sys.exit(1)
        tid, srd, bseq, mseq = ary
        
        dna = Seq(bseq)
        if srd == '-': dna = dna.reverse_complement()
        n1 = dna.translate().count("*")
        bseq = dna.translate(to_stop = True)
        bpro = str(bseq)
        
        dna = Seq(mseq)
        if srd == '-': dna = dna.reverse_complement()
        n2 = dna.translate().count("*")
        mseq = dna.translate(to_stop = True)
        mpro = str(mseq)
        fho.write("%s\t%d\t%d\t%s\t%s\n" % (tid, n1, n2, bpro, mpro))
        rcd = SeqRecord(mseq, id = tid, description = '')
        rcds.append(rcd)
    SeqIO.write(rcds, fhs, "fasta")

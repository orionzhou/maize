#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

def read_seqid_map(fmap):
    fhi = open(fmap, "r")
    mdict = dict()
    for line in fhi:
        line = line.strip("\n")
        ochr, length, nchr = line.split("\t")
        if ochr.count("~") > 0:
            ochr = ochr.split("~")[0]
        mdict[ochr] = nchr
    fhi.close()
    return mdict
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'gff'
    )
    parser.add_argument(
        'fi', help = 'input gff'
    )
    parser.add_argument(
        'fo', help = 'output gff'
    )
    parser.add_argument(
        '--map', dest = "fmap", help = 'seqid mapping file'
    )
    args = parser.parse_args()
    
    mdict = read_seqid_map(args.fmap)
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        line = line.strip("\n")
        if line.startswith("#"):
            print >>fho, line
            continue
        ary = line.split("\t")
        if len(ary) < 9:
            continue
        seqid, src, type, beg, end, score, srd, phase, desc = ary
#        if type in ["contig", "match", "match_part", "protein_match", "expressed_sequence_match"]:
        if type not in ["gene", "mRNA", "CDS", 'five_prime_UTR', 'three_prime_UTR']:
            continue
        if not mdict.has_key(seqid):
            print "%s not in seqid mapping file" % seqid
            sys.exit(1)
        seqid = mdict[seqid]
        ary = [seqid, src, type, beg, end, score, srd, phase, desc]
        print >>fho, "\t".join(ary) 
    fhi.close()
    fho.close()

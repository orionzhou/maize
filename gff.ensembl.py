#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

def parse_desc(desc):
    ary = desc.split(";")
    keys, vals, rdic = list(), list(), dict()
    for ele in ary:
        key, val = ele.split('=')
        ary2 = val.split(":")
        if len(ary2) > 1:
            val = ary2[1]
        keys.append(key)
        vals.append(val)
        rdic[key] = val
    return keys, vals, rdic
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
    args = parser.parse_args()
    
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
        seqid, src, featype, beg, end, score, srd, phase, desc = ary
        keys, vals, rdic = parse_desc(desc)
#        if type in ["contig", "match", "match_part", "protein_match", "expressed_sequence_match"]:
        if not seqid in ["%d" % x for x in range(1, 11)]:
            continue
        if featype == "transcript":
            if rdic.has_key('biotype') and rdic['biotype'] == 'protein_coding':
                featype = "mRNA"
            else:
                continue
        if featype not in ["gene", "mRNA", "CDS", 'five_prime_UTR', 'three_prime_UTR']:
            continue
        desc = ";".join(["=".join(x) for x in zip(keys, vals)])
        ary = [seqid, src, featype, beg, end, score, srd, phase, desc]
        print >>fho, "\t".join(ary) 
    fhi.close()
    fho.close()

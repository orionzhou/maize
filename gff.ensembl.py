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
    seqids1 = ["%d" % x for x in range(1, 11)]
    for line in fhi:
        line = line.strip("\n")
        if line.startswith("#"):
            fho.write(line + "\n")
            continue
        ary = line.split("\t")
        if len(ary) < 9:
            continue
        seqid, src, featype, beg, end, score, srd, phase, desc = ary
        keys, vals, rdic = parse_desc(desc)
        if seqid not in seqids1 and not seqid.startswith("chr"):
            continue
        if featype == "transcript":
            if 'biotype' in rdic and rdic['biotype'] == 'protein_coding':
                featype = "mRNA"
            else:
                continue
        if featype not in ["gene", "mRNA", "CDS", 'five_prime_UTR', 'three_prime_UTR']:
            continue
        desc = ";".join(["=".join(x) for x in zip(keys, vals)])
        ary = [seqid, src, featype, beg, end, score, srd, phase, desc]
        fho.write("\t".join(ary) + "\n")
    fhi.close()
    fho.close()

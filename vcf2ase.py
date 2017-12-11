#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert vcf to tsv'
    )
    parser.add_argument(
        'fi', help = 'input file (.vcf)'
    )
    parser.add_argument(
        'fo', help = 'output file (.tsv)'
    )
    args = parser.parse_args()

    if args.fi and args.fi != '-' and args.fi != 'stdout':
        sys.stdin = open(args.fi, "r")
    fho = open(args.fo, "w")
    lsth = ["chr", "pos", "ref", "alt", "qual", "depth", "dpr", "dpa"]
    fho.write("\t".join(lsth) + "\n")

    for line in sys.stdin:
        if line.startswith("#"): continue
        ps = line.strip().split("\t")
        chrom, pos, vid, ref, alt, qual, flt, info, fmt, call = ps
        tdic = dict()
        for tagstr in info.split(";"):
            if "=" not in tagstr: continue
            tag, val = tagstr.split("=")
            tdic[tag] = val

        dp, dpr, dra = 0, 0, 0
        if 'DP' in tdic:
            dp = int(tdic['DP'])
        if 'DP4' in tdic:
            dp4 = [int(x) for x in tdic['DP4'].split(",")]
            dpr = dp4[0] + dp4[1]
            dpa = dp4[2] + dp4[3]
        #if dp != dpr + dpa:
        #    print("%s:%s %d <> %d + %d" % (chrom, pos, dp, dpr, dpa))
        lst = [chrom, pos, ref, alt, qual, dp, dpr, dpa]
        fho.write("\t".join([str(x) for x in lst]) + "\n")
    fho.close()
    

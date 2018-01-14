#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert inbred vcf to hybrid vcf'
    )
    parser.add_argument(
        'fi', help = 'input file (.vcf)'
    )
    parser.add_argument(
        'fo', help = 'output file (.vcf)'
    )
    args = parser.parse_args()

    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith("##fileformat") or \
                line.startswith("##reference") or \
                line.startswith("##FORMAT=<ID=GT") or \
                line.startswith("##contig"):
            fho.write(line)
            continue
        elif line.startswith("##"):
            continue
        row = line.strip("\n").split("\t")
        if row[0] == "#CHROM":
            row[9] = "B73xMo17"
        else:
            row[6] = "PASS"
            row[7] = "."
            row[8] = "GT"
            row[9] = "0|1"
        fho.write("\t".join(row) + "\n")
    fhi.close()
    fho.close()


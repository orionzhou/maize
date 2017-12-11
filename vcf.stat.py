#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'calculate stats for each site'
    )
    parser.add_argument(
        'fi', help = 'input file (VCF)'
    )
    parser.add_argument(
        'fo', nargs = '?', help = 'output file (tsv)'
    )
    args = parser.parse_args()

    if args.fo and args.fo != '-' and args.fo != 'stdout':
        sys.stdout = open(args.fo, "w")
    print "\t".join(["chr", "pos", "nalt", "rsize", "asize", "nsam", "aaf", "nucdiv"])

    vcf_reader = vcf.Reader(open(args.fi, 'r'))
    for rcd in vcf_reader:
        num_chroms = float(2.0 * rcd.num_called)
        nucl_diversity = float(num_chroms / (num_chroms - 1.0)) * rcd.heterozygosity
        print "\t".join(map(str, [rcd.CHROM, rcd.POS,
            len(rcd.ALT), len(rcd.REF), len(rcd.ALT[0]), 
            rcd.num_called, rcd.aaf[0], nucl_diversity]))
    

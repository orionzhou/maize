#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert VCF to Tbl file'
    )
    parser.add_argument(
        'fi', help = 'input file (VCF)'
    )
    parser.add_argument(
        'fo', nargs = '?', help = 'output file (Tbl)'
    )
    args = parser.parse_args()

    if args.fo and args.fo != '-' and args.fo != 'stdout':
        sys.stdout = open(args.fo, "w")

    vcf_reader = vcf.Reader(open(args.fi, 'r'))
    print "\t".join(['chr', 'pos', 'alt', 'gt', 'rd', 'qual', 'mapqual'])
    for rcd in vcf_reader:
        samples = rcd.samples
        assert len(samples) == 1, "%d samples in %s" % (len(samples), rcd)
        sample = samples[0]
        alt = rcd.ALT[0]
        print "\t".join(map(str, [rcd.CHROM, rcd.POS, alt,
            sample.gt_type, sample.data.DP, rcd.QUAL, rcd.INFO['MQ']]))
    

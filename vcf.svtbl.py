#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert SV VCF to Tbl'
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
    print "\t".join(['chr', 'beg', 'end'])
    for rcd in vcf_reader:
        samples = rcd.samples
        alts = rcd.ALT
        assert len(samples) == 1, "%d samples in %s" % (len(samples), rcd)
        assert len(alts) == 1, "%d alts in %s" % (len(alts), rcd)
        sample = samples[0]
        alt = rcd.ALT[0]
        beg = rcd.POS
        rlen = len(rcd.REF)
        alen = len(alt)
        if rlen == 1:
            continue
        print "\t".join(map(str, [rcd.CHROM, beg + 1, beg + rlen - 1]))
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def main(args):
    pre = args.fo
    sh(f"rm {pre}.*")
    cmd = ''
    if args.s2 == '':
        samples = args.s1
        cmd = (
            f'bcftools view -s {samples} -Ou {args.vcf} |',
            f'bcftools filter -i \'N_PASS(GQ>={args.gq} & GT!="mis")=1 && FORMAT/GT[0]="AA"\' |'
            'bcftools annotate -x QUAL,INFO,^FORMAT/GT |'
            'bioawk -t \'{if(!/^#/){$10=substr($10,1,1)"|"substr($10,1,1)}; print}\' |',
            f'bcftools reheader -h {args.header} > {pre}.vcf'
        )
    else:
        samples = f"{args.s1},{args.s2}"
        ref_rule = ''
        if args.s1.endswith("B73"):
            ref_rule = '&& FORMAT/GT[0]="RR"'
        elif args.s2.endswith("B73"):
            ref_rule = '&& FORMAT/GT[1]="RR"'
        cmd = (
            f'bcftools view -s {samples} -Ou {args.vcf} |',
            f'bcftools filter -i \'N_PASS(GQ>={args.gq} & GT!="mis")=2 && N_PASS(GT="RR")=1 && N_PASS(GT="AA")=1 {ref_rule}\' |',
            'bcftools annotate -x QUAL,INFO,^FORMAT/GT |',
            'bioawk -t \'{if(!/^#/){$10=substr($10,1,1)"|"substr($11,1,1)}; $11=""; print}\' |',
            f'bcftools reheader -h {args.header} > {pre}.vcf'
        )
    sh(''.join(cmd))
    sh(f"bgzip {pre}.vcf")
    sh(f"bcftools index -t {pre}.vcf.gz")
    sh(f"bcftools view {pre}.vcf.gz -Ob -o {pre}.bcf")
    sh(f"bcftools index {pre}.bcf")
    sh(f"bcftools stats -s - {pre}.bcf > {pre}.txt")
    sh(f"bcftools query -f '%CHROM\\t%POS\t[%TGT]\\n' {pre}.bcf | sed 's/|/,/' > {pre}.tsv")

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "create variant files for nextflow/ASE pipeline"
    )

    ps.add_argument('s1', help = 'sample 1')
    ps.add_argument('fo', help = 'output file prefix')
    ps.add_argument('--s2', default='', help = 'sample 2')
    ps.add_argument('--vcf', default=f"{os.environ['s3']}/zhoup-nfo/zm.vt10/04.snp.vcf.gz", help='input (joint) VCF')
    ps.add_argument('--gq', default=30, help='minimum Genotype Quality score')
    ps.add_argument('--header', default=f"{os.environ['maize']}/assets/dummy_header.vcf", help='vcf header')

    args = ps.parse_args()
    main(args)


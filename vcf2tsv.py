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

    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    vcf_reader = vcf.Reader(fhi)
    lst1 = ["DP", "QD", "FS", "MQ", "MQRankSum", "ReadPosRankSum", "SOR"]
    lsth = ['chr', 'pos', 'ref', 'alt', 'IS_SNP', 'PASS', 'QUAL'] + lst1
    lst2 = ["AD", "DP", "GQ"]
    lstb = ['GT'] + lst2
    head = 1
    for rcd in vcf_reader:
        n_sm = len(rcd.samples)
        if head == 1:
            head = 0
            if n_sm == 1:
                fho.write("\t".join(lsth + lstb) + "\n")
            else:
                sm_names = [x.sample for x in rcd.samples]
                lstbe = []
                for sm_name in sm_names:
                    lstbe1 = [i+j for i,j in zip([sm_name]*len(lstb), lstb)]
                    lstbe += lstbe1
                fho.write("\t".join(lsth + lstbe) + "\n")
        alts = ",".join(map(str, rcd.ALT))
        alt = rcd.ALT[0]
        filt = rcd.FILTER
        flagpass = 0
        if filt is None or len(filt) == 0:
            flagpass = 1
        val1, val2 = [], []
        for k in lst1:
            v = ''
            if k in rcd.INFO:
                v = rcd.INFO[k]
            val1.append(v)
        valh = [rcd.CHROM, rcd.POS, rcd.REF, alts, int(rcd.is_snp), 
                flagpass, rcd.QUAL] + val1
        sms = rcd.samples
        valb = []
        for sm in sms: # need to address multiple alleles
            val2 = [getattr(sm.data, k, '') for k in lst2]
            if val2[0] is not None and val2[0] != '':
                val2[0] = val2[0][1]
            valb1 = [sm.gt_type] + val2
            valb += ['' if x is None else str(x) for x in valb1] 
        fho.write("\t".join(map(str, valh + valb)) + "\n")
    

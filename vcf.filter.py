#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import vcf

def read_sites(fi):
    sites = set()
    fhi = open(fi, "r")
    for line in fhi:
        seqid, pos = line.strip().split("\t")
        locus = "%s_%s" % (seqid, pos)
        sites.add(locus)
    return sites
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'filter vcf'
    )
    parser.add_argument(
            'fi', help = 'input file (.vcf)'
    )
    parser.add_argument(
            'fo', help = 'output file (.vcf)'
    )
    parser.add_argument(
            'fx', help = 'sites to exclude (.tsv)'
    )
    args = parser.parse_args()
    sites = read_sites(args.fx)

    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    vcfr = vcf.Reader(fhi)
    #vcfw = vcf.Writer(fho, vcfr)
    for rcd in vcfr:
        n_sm = len(rcd.samples)
        sms = rcd.samples
        locus = "%s_%s" % (rcd.CHROM, rcd.POS)
        if locus in sites:
            continue
        sm = sms[0]
        gt = sm.gt_type
        qd = rcd.INFO['QD'] if 'QD' in rcd.INFO else None
        fs = rcd.INFO['FS'] if 'FS' in rcd.INFO else None
        mq = rcd.INFO['MQ'] if 'MQ' in rcd.INFO else None
        mqrs = rcd.INFO['MQRankSum'] if 'MQRankSum' in rcd.INFO else None
        rprs = rcd.INFO['ReadPosRankSum'] if 'ReadPosRankSum' in rcd.INFO else None
        sor = rcd.INFO['SOR'] if 'SOR' in rcd.INFO else None
        flagpass = (gt is None or gt == 2) and (
                (rcd.is_snp &
                    (qd is None or qd >= 2) &
                    (fs is None or fs <= 60) &
                    (mq is None or mq >= 40) &
                    (mqrs is None or mqrs >= -12.5) &
                    (rprs is None or rprs >= -8) & 
                    (sor is None or sor <= 4)
                ) or 
                (rcd.is_indel & 
                    (qd is None or qd >= 2) &
                    (fs is None or fs <= 200) &
                    (rprs is None or rprs >= -20) &
                    (sor is None or sor <= 10)
                )
        )
        alts = ",".join(map(str, rcd.ALT))
        alt = str(rcd.ALT[0])
        if flagpass:
            if rcd.is_snp:
                fho.write("%s\t%s\t%s\t%d\t%s\n" % (locus, 'single', rcd.CHROM, rcd.POS-1, alt))
            elif rcd.is_deletion:
                fho.write("%s\t%s\t%s\t%d\t%d\n" % (locus, 'deletion', rcd.CHROM, rcd.POS-1, len(rcd.REF)-1))
            else:
                fho.write("%s\t%s\t%s\t%d\t%s\n" % (locus, 'insertion', rcd.CHROM, rcd.POS-1, alt[1:]))
            #vcfw.write_record(rcd)
    #vcfw.close()
    

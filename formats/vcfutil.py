#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import vcf

from maize.apps.base import sh, mkdir
from maize.formats.base import must_open

dna_comp_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def revcomp(dna):
    return ''.join([dna_comp_dict[base] for base in dna[::-1]])

def stat(args):
    fhi = must_open(args.fi)
    print("\t".join(["chr", "pos", "nalt", "rsize", "asize", "nsam", "aaf", "nucdiv"]))

    vcf_reader = vcf.Reader(fhi)
    for rcd in vcf_reader:
        num_chroms = float(2.0 * rcd.num_called)
        nucl_diversity = float(num_chroms / (num_chroms - 1.0)) * rcd.heterozygosity
        print("\t".join(map(str, [rcd.CHROM, rcd.POS, \
            len(rcd.ALT), len(rcd.REF), len(rcd.ALT[0]), \
            rcd.num_called, rcd.aaf[0], nucl_diversity])))

def vcf_filter(args):
    sites = set()
    if args.exclude:
        for line in must_open(args.exclude):
            seqid, pos = line.strip().split("\t")
            locus = "%s_%s" % (seqid, pos)
            sites.add(locus)
    
    vcfr = vcf.Reader(must_open(args.fi))
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
            #vcfw.write_record(rcd)
            if rcd.is_snp:
                print("%s\t%s\t%s\t%d\t%s" % (locus, 'single', rcd.CHROM, rcd.POS-1, alt))
            elif rcd.is_deletion:
                print("%s\t%s\t%s\t%d\t%d" % (locus, 'deletion', rcd.CHROM, rcd.POS-1, len(rcd.REF)-1))
            else:
                print("%s\t%s\t%s\t%d\t%s" % (locus, 'insertion', rcd.CHROM, rcd.POS-1, alt[1:]))

def vcf2fas(args):
    from Bio import SeqIO

    fhi = must_open(args.fi)
    seqdic = SeqIO.index(args.fs, "fasta")
    vcfr = vcf.Reader(filename = args.fv, compressed = True)
     
    for line in fhi:
        line = line.strip("\n")
        ps = line.split("\t")
        if len(ps) < 4:
            print("<4 fields: %s" % line)
            continue
        seqid, beg, end, note = ps[0:4]
        beg, end = int(beg), int(end)
        chrstr = seqdic[seqid].seq
        seqstr = chrstr[beg-1:end]
        rcds = vcfr.fetch(seqid, start = beg-1, end = end)
        for rcd in rcds:
            alts = ",".join(map(str, rcd.ALT))
            alt = rcd.ALT[0]
            #logging.debug("%s\t%s\t%s\t%s" % (rcd.CHROM, rcd.POS, rcd.REF, alts))
        #logging.debug("%s\t%s\t%s\t%s" % (seqid, beg, end, note))
        print("%s\t%s\t%s\t%s\t%s" % (seqid, beg, end, note, seqstr))

def vcf2bed(args):
    # bcftools query -f '%CHROM\t%POS0\t%END\t%TYPE:%ALT\t0|1\n' Mo17.vcf -o x.bed
    fhi = must_open(args.fi)
    for line in fhi:
        if line.startswith("#"):
            continue
        row = line.strip("\n").split("\t")
        sid, pos, ref, alt, phase = row[0], row[1], row[3], row[4], row[9]
        pos = int(pos)
        vnt = "M:%s" % alt
        if len(ref) > 1:
            assert len(alt) == 1, "error: %s" % line
            vnt = "D:%d" % (len(ref) - 1)
            pos = pos + 1
            if args.noindel:
                continue
        elif len(alt) > 1:
            assert len(ref) == 1, "error: %s" % line
            vnt = "I:%s" % alt[1:]
            if args.noindel:
                continue
        else:
            assert len(ref) == 1 and len(alt) == 1, "error: %s" % line
        print("\t".join([sid, str(pos-1), str(pos), vnt, phase]))

def vcf2tsv(args):
    vcf_reader = vcf.Reader(fsock=must_open(args.fi))
    #vcf_reader.fetch("B02", 50001, 100000)
    #print(str(vcf_reader))
    #sys.exit()
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
                print("\t".join(lsth + lstb))
            else:
                sm_names = [x.sample for x in rcd.samples]
                lstbe = []
                for sm_name in sm_names:
                    lstbe1 = [i+'_'+j for i,j in zip([sm_name]*len(lstb), lstb)]
                    lstbe += lstbe1
                print("\t".join(lsth + lstbe))
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
        print("\t".join(map(str, valh + valb)))

def ase(args):
    lsth = ["chr", "pos", "ref", "alt", "qual", "depth", "dpr", "dpa"]
    print("\t".join(lsth))
    for line in must_open(args.fi):
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
        print("\t".join([str(x) for x in lst])) 

def hybrid(args):
    for line in must_open(args.fi):
        line = line.strip()
        if line.startswith("##fileformat") or \
                line.startswith("##reference") or \
                line.startswith("##FORMAT=<ID=GT") or \
                line.startswith("##contig"):
            print(line)
            continue
        elif line.startswith("##"):
            continue
        row = line.split("\t")
        if row[0] == "#CHROM":
            row[9] = "B73xMo17"
        else:
            row[6] = "PASS"
            row[7] = "."
            row[8] = "GT"
            row[9] = "0|1"
        print("\t".join(row))
    fhi.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'vcf utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')
    
    sp1 = sp.add_parser("stat", help = "report stats for each variant")
    sp1.add_argument('fi', help = 'input vcf file')
    sp1.set_defaults(func = stat)

    sp1 = sp.add_parser("filter", help = "filter sites")
    sp1.add_argument('fi', help = 'input vcf file')
    sp1.add_argument('--exclude', default = None, help = 'sites to exclude (.tsv)')
    sp1.set_defaults(func = vcf_filter)
    
    sp1 = sp.add_parser("2fas", help = "make an alternate fasta using vcf")
    sp1.add_argument('fi', help = 'interval file (.tsv)')
    sp1.add_argument('fv', help = 'vcf file (.vcf)')
    sp1.add_argument('fs', help = 'reference sequence file (.fas)')
    sp1.set_defaults(func = vcf2fas)
 
    sp1 = sp.add_parser("2tsv", help = "vcf -> tsv")
    sp1.add_argument('fi', help = 'input vcf file')
    sp1.set_defaults(func = vcf2tsv)

    sp1 = sp.add_parser("ase", help = "report allele proportion (ASE) for each variant")
    sp1.add_argument('fi', help = 'input vcf file')
    sp1.set_defaults(func = ase)

    sp1 = sp.add_parser("hybrid", help = "convert inbred vcf to hybrid vcf")
    sp1.add_argument('fi', help = 'input vcf file')
    sp1.set_defaults(func = hybrid)

    sp1 = sp.add_parser("2bed", help = "extract variant locations")
    sp1.add_argument('fi', help = 'input vcf file')
    sp1.add_argument('--noindel', action = 'store_true', help = 'remove InDel')
    sp1.set_defaults(func = vcf2bed)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

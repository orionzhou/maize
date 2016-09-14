#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import argparse
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
def read_cfg(fc):
    (orgs, fis) = ([], [])
    fhc = open(fc, "r")
    for line in fhc:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        (org, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            print os.access(fi, os.F_OK)
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
    fhc.close()
    return (orgs, fis)
def merge_seqs(fis, fids, fo):
    print "  merging input files to %s" % fo
    seqs = []
    for i in range(0,len(fids)):
        handle = 0
        if (fis[i].endswith(".gz")):
            handle = gzip.open(fis[i], "rb")
        else:
            handle = open(fis[i], "rU")
        seq_it = SeqIO.parse(handle, "fasta")
        handle.close

        seqs1 = [SeqRecord(rcd.seq, id = fids[i] + "|" + rcd.id,
            description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        'org', nargs = "?", help = 'accession'
    )
    parser.add_argument(
        'outfile', nargs = "?", help = 'output fasta'
    )
    args = parser.parse_args()
    #(fc, fo) = (args.cfgfile, args.outfile)
    
    org = "HM101"
    dirw = op.join(os.environ['misc2'], "nrgenome", org)
    os.chdir(dirw)
    
    #os.system("cat 11.blat/*.psl > 12.psl")
    
    fho = open("14.tbl", 'w')
    for res in SearchIO.parse('12.psl', 'blat-psl'):
        print res.id
        for hit in res:
            print "\t%s" % hit.id
            for hsp in hit:
                print "\t\t[%d-%d] [%d-%d] %d %g" % (hsp.query_start, hsp.query_end, hsp.hit_start, hsp.hit_end, hsp.ident_num, hsp.ident_pct)
    fho.close()

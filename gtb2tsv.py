#!/usr/bin/env python
import os
import os.path as op
import sys
from CoordSys import *

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'gff to tsv'
    )
    parser.add_argument(
            'fi', help = 'input gff'
    )
    parser.add_argument(
            'fo', help = 'output tsv'
    )
    args = parser.parse_args()
    
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    fho.write("gid\ttid\tttype\tetype\tseqid\tbeg\tend\tsrd\n")
    for line in fhi:
        line = line.strip("\n")
        if line.startswith("#") or line.startswith("id"):
            continue
        ary = line.split("\t")
        if len(ary) < 18:
            print("less than 18 columns:\n%s" % line)
            continue
        tid, gid, seqid, tbeg, tend, srd, \
                locES, locIS, locCS, loc5S, loc3S, phase, \
                src, conf, cat1, cat2, cat3, note = ary
        tbeg, tend = int(tbeg), int(tend)
        if cat1 == 'mRNA':
            assert locCS, "no CDS for %d" % tid
        else:
            assert locES, "no exon for %d" % tid
        ldic = { 'exon': locES, 'cds': locCS, \
                'utr5': loc5S, 'utr3': loc3S, 'intron':locIS }
        for etype, locS in ldic.items():
            if not locS:
                continue
            for rbeg, rend in locStr2Ary(locS):
                beg, end = 0, 0
                if srd == "-":
                    beg, end = tend - rend + 1, tend - rbeg + 1
                else:
                    assert srd == '+', "unknown strand: %s for %s" % (srd, tid)
                    beg, end = tbeg + rbeg - 1, tbeg + rend - 1
                fields = [gid, tid, cat1, etype, seqid, str(beg), str(end), srd]
                fho.write("\t".join(fields)+"\n") 
    fhi.close()
    fho.close()

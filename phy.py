#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'make multiple sequence alignment and build ML tree'
    )
    parser.add_argument(
        'fi', help = 'input sequence file (fasta)'
    )
    parser.add_argument(
        'fo', help = 'output prefix'
    )
    parser.add_argument(
        '--seqtype', default = 'nt', help = 'type of sequence [nt/aa]'
    )
    args = parser.parse_args()
    (fi, fo, seqtype) = (args.fi, args.fo, args.seqtype)

    os.system("seq.desc.py %s %s.tsv" % (fi, fi))
    nproc = int(os.environ['nproc'])
    if seqtype == 'aa':
        cmd = "clustalo -i %s -o %s.aln --outfmt=clu --force --use-kimura --threads=%d" % (fi, fo, nproc)
    elif seqtype == 'nt':
        #cmd = "clustalw2 -INFILE=%s -ALIGN -TYPE=DNA -OUTORDER=INPUT -OUTFILE=%s.aln" % (fi, fo)
        cmd = "muscle -in %s -out %s.aln -clw" % (fi, fo)
    else:
        print "unknown seqtype: %s" % seqtype
        sys.exit(1)
    print cmd
    os.system(cmd)
    ##os.system("trimal -in %s.aln -out %s.trimmed.aln -automated1" % (fo, fo))
    os.system("aln.conv.py --fmt 2 %s.aln %s.phy" % (fo, fo))
    #os.system("phyml -i %s.phy -d %s -o tlr --quiet" % (fo, seqtype))
    #os.system("mv %s.phy_phyml_tree.txt %s.nwk" % (fo, fo))
    #os.system("rm %s.phy %s.phy_phyml*" % (fo, fo))
 

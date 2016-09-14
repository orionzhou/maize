#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'run PhyML on input alignment (.aln)'
    )
    parser.add_argument(
        'fi', help = 'input alignment (.aln)'
    )
    parser.add_argument(
        'fo', help = 'output prefix'
    )
    parser.add_argument(
        '--seqtype', default = 'nt', help = 'type of sequence [nt/aa]'
    )
    args = parser.parse_args()
    (fi, fo, seqtype) = (args.fi, args.fo, args.seqtype)

    nproc = int(os.environ['nproc'])
    os.system("aln.conv.py --fmt 2 %s %s.phy" % (fi, fo))
    os.system("phyml -i %s.phy -d %s -o tlr --quiet" % (fo, seqtype))
    os.system("mv %s.phy_phyml_tree.txt %s.nwk" % (fo, fo))
    os.system("rm %s.phy %s.phy_phyml*" % (fo, fo))
 

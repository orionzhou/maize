#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'unknown'
    )
    parser.add_argument(
        'org', nargs = "?", default = "HM034", help = 'genome to valid'
    )
    args = parser.parse_args()
    org = args.org

    dirw = op.join(os.environ['misc3'], "alpaca/34.sfams")
    os.chdir(dirw)

    sfams = ['CRP0355', 'CRP3710', 'CRP4180']
    for sfam in sfams:
        for opt in ['alps', 'alpc']:
            pre = "%s.%s" % (sfam, opt)
            os.system("seqret.pl -d ../32.pro.fas -b %s.bed -o %s.fas" % (pre, pre))
            os.system("clustalo -i %s.fas -o %s.aln --outfmt=clu --force --use-kimura --threads=4" % (pre, pre))
            os.system("aln2phy.pl -i %s.aln -o %s.phy -l 30" % (pre, pre))
            os.system("phyml -i %s.phy -d aa" % pre)
            os.system("mv %s.phy_phyml_tree.txt %s.phy.nwk" % (pre, pre))
            os.system("rm %s.phy_phyml*" % pre)
            #sys.exit()
    

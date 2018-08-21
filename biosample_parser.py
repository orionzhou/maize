#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
from bs4 import BeautifulSoup

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fasta utilities'
    )
    parser.add_argument('fi', help = 'input biosample xml')
    args = parser.parse_args()

    fhi = open(args.fi, 'r')
    soup = BeautifulSoup(fhi, 'xml')
    print("BioSample\tTitle\ttissue\tgenotype")
    for bs in soup.BioSampleSet.children:
        if len(bs) == 1: continue
        bsid = bs['accession']
        title = bs.Description.Title.string.strip("\"")
        tissue = ''
        genotype = ''
        if bs.Attributes:
            for x in bs.Attributes.children:
                if x.string.strip() == '': continue
                if x['attribute_name'] == 'tissue':
                    tissue = x.string.strip()
                elif x['attribute_name'] == 'genotype':
                    genotype = x.string.strip()
        print("\t".join([bsid, title, tissue, genotype]))

 

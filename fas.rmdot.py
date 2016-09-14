#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
from string import maketrans

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'replace periods (.) in an fasta file by dash (-)'
    )
    parser.add_argument(
        'fi', help = 'input file (.fas)'
    )
    parser.add_argument(
        'fo', help = 'output file (.fas)'
    )
    args = parser.parse_args()
    
    tt = maketrans(".", "-")
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith('>'):
            fho.write(line)
        else:
            fho.write(line.translate(tt))
    fhi.close()
    fho.close()

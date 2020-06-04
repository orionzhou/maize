#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import os.path as op
from Bio import motifs

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'convert cis_bp PWM file to chen format'
    )
    parser.add_argument('fi', help = 'input cis_bp PWM file')
    parser.add_argument('--id', default = '', help = 'motif id [file basename]')
    args = parser.parse_args()

    mid = args.id
    if mid == '':
        mid = op.basename(op.splitext(args.fi)[0])

    fhi = open(args.fi, 'r')
    for m in motifs.parse(fhi, 'minimal'):
        print("%s\t%s" % (m.name, m.evalue))
    fhi.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path as op
import argparse
sys.path.append("/home/springer/zhoux379/source/genie3/GENIE3_python")
from GENIE3 import *
from numpy import loadtxt

if __name__ == '__main__':
    dirw = '/home/springer/zhoux379/data/misc2/grn23/62.genie3'
    fi = op.join(dirw, "01.matrix.tsv")
    fr = op.join(dirw, "11.TF.txt")
    fo = op.join(dirw, "31.tsv")
    
    data = loadtxt(fi, skiprows=1)
    fhi = open(fi, "r")
    gids = fhi.readline()
    fhi.close()
    gids = gids.rstrip('\n').split('\t')

    rids = [line.strip() for line in open(fr, 'r')]
    VIM = genie3(data, gene_names = gids, regulators = rids)
    get_link_list(VIM, gene_names = gids, regulators = rids, file_name = fo)

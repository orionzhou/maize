#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import numpy as np
import argparse
import networkx as nx

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'extract all connected components of a graph'
    )
    parser.add_argument(
        'fi', help = 'input file (tbl)'
    )
    parser.add_argument(
        'fo', help = 'output file (tbl)'
    )
    args = parser.parse_args()
    
    fho = open(args.fo, "w")
    
    ary = np.genfromtxt(args.fi, names = True, dtype = None)
    edges = map(lambda ary: (ary[0], ary[1]), ary)
    G = nx.Graph()
    G.add_edges_from(edges)
    conn = nx.connected_components(G)
    
    fho.write("clu\tid\n")
    cluster = 1
    for names in conn:
        for name in names:
            fho.write("%d\t%s\n" % (cluster, str(name, 'utf-8')))
        cluster += 1


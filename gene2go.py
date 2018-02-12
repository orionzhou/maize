#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
from goatools import obo_parser
from download import download_file

def get_lev2_ancestors(G, node):
    ancs = nx.descendants(G, node)
    ancs_lev1 = set(x for x in ancs if G.out_degree(x) == 0)
    ancs_lev2 = set()
    for x in ancs:
        for y in G.successors(x):
            if y in ancs_lev1:
                ancs_lev2.add(x)
                break
    return ancs_lev2
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Generate Gene2Go file'
    )
    parser.add_argument(
            'fi', help = 'input file'
    )
    parser.add_argument(
            'fg', help = 'output gene2go table'
    )
    parser.add_argument(
            'fd', help = 'output GO discription table'
    )
    parser.add_argument(
            '--no_propagate', default = False, action = 'store_true',
            help = 'Do not propogate counts to parent terms'
    )
    parser.add_argument(
            '--download', default = False, action = 'store_true',
            help = 'Download/update go-basic.obo'
    )
    parser.add_argument(
            '--obo', default = '/home/springer/zhoux379/data/db/go/go-basic.obo',
            help = 'Path of go-basic.obo'
    )
    args = parser.parse_args()
    (fi, fg, fd) = (args.fi, args.fg, args.fd)

    if args.download:
        if op.isfile(args.obo):
            os.remove(args.obo)
        url = "http://geneontology.org/ontology/go-basic.obo"
        download_file(url, args.obo)
    assert op.isfile(args.obo), "go-basic.obo not found: %s" % args.obo

    dag = obo_parser.GODag(args.obo)
    fhd = open(fd, "w")
    fhd.write("goid\tnamespace\tlevel\tdepth\tname\n")
    for go, node in dag.items():
        if node.is_obsolete:
            continue
        fhd.write("%s\t%s\t%d\t%d\t%s\n" % (go, node.namespace, 
            node.level, node.depth, node.name))
    fhd.close()
    
    fhi = open(fi, "r")
    fhg = open(fg, "w")
    for line in fhi:
        line = line.strip("\n")
        gid, gostr = line.split("\t")
        gos = gostr.split(";")
        gos = [x for x in gos if x]
        gset = set()
        for go in gos:
            if not go in dag:
                continue
            ancs = dag[go].get_all_parents()
            if len(ancs) == 0:
                continue
            gset.add(go)
            if not args.no_propagate:
                gset |= ancs
        if len(gset) > 0:
            for go in gset:
                fhg.write("%s\t%s\n" % (gid, go))
    fhi.close()
    fhg.close()

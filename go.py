#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import argparse
import networkx as nx
import obonet

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
    parser = argparse.ArgumentParser(
        description = 'make multiple sequence alignment and build ML tree'
    )
    parser.add_argument(
        'fi', help = 'input file'
    )
    parser.add_argument(
        'fo', help = 'output file'
    )
    args = parser.parse_args()
    (fi, fo) = (args.fi, args.fo)

    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    url = "http://geneontology.org/ontology/go-basic.obo"
    G = obonet.read_obo(url)
    print(len(G))
    id_to_name = {id_: data['name'] for id_, data in G.nodes(data=True)}
    
    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        line = line.strip("\n")
        gid, gostr = line.split("\t")
        gos = gostr.split(";")
        gos = [x for x in gos if x]
        gos2 = set()
        for go in gos:
            if not G.has_node(go):
                continue
            ancs_lev2 = get_lev2_ancestors(G, go)
            gos2 = gos2 | ancs_lev2
        gos2 = gos
        if len(gos2) > 0:
            for go in gos2:
                if go in id_to_name:
                    fho.write("%s\t%s\t%s\n" % (gid, go, id_to_name[go]))
    fhi.close()
    fho.close()

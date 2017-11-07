#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path as op
import networkx as nx
#from time import clock, time

def trim_pathname(name):
    name = name.strip("\"").replace("&","").replace(";","")
    name = name.replace("<i>","").replace("</i>","")
    name = name.replace("<sub>","").replace("</sub>","")
    name = name.replace("<sup>","").replace("</sup>","")
    return name
def get_parents(G, node, palist):
    pas = G.predecessors(node)
    if len(pas) > 0:
        pa = pas[0]
        palist.append(pa)
        get_parents(G, pa, palist)
def corncyc2tsv():
    dirw = '/home/springer/zhoux379/data/genome/Zmays_v4/corncyc'
    fi = op.join(dirw, "query-results.tsv")
    fo = op.join(dirw, "01.tsv")
    fhi = open(fi, "r")
    fho = open(fo, "w")
    pdic = dict()
    nodes, edges = set(), []
    for line in fhi:
        ps = line.strip().split("\t")
        if len(ps) < 2:
            continue
        else:
            name, gidstr, pastr = ps[0:3]
            name = trim_pathname(name)
            gids = gidstr.strip("()").split(" ")
            gids = [x.strip("\"").split("_")[0] for x in gids]
            gids = [x for x in gids if x != '']
            pas = pastr.strip("()").split("\" \"")
            pas = [trim_pathname(x.strip("\"")) for x in pas]
            if pas[0] == "Pathways": continue
            if name in nodes:
                print("appeared >1 times: %s" % name)
                continue
            pdic[name] = gids
            nodes.add(name)
            if len(pas) > 0:
                edges.append((pas[0], name))
#            fho.write("%s\t%s\t%s\n" % (name, "|".join(pas), ",".join(gids)))
#            if len(gids) > 0:
#                for gid in gids:
#                    fho.write("%s\t%s\n" % (name, gid))
    print(len(nodes))
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    #assert nx.is_directed_acyclic_graph(G), "not a DAG"
    allnodes = G.nodes()
    leaves = []
    for node in allnodes:
        if len(G.successors(node)) == 0:
            gids = pdic[node]
            if len(gids) == 0: continue
            ancs = []
            get_parents(G, node, ancs)
            ancs = list(reversed(ancs))
            (anc1, anc2, anc3, anc4) = (node, node, node, node)
            if len(ancs) > 0:
                anc1 = ancs[0]
            if len(ancs) > 1:
                anc2 = ancs[1]
            if len(ancs) > 2:
                anc3 = ancs[2]
            if len(ancs) > 3:
                anc4 = ancs[3]
            fho.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (node, anc1, anc2, anc3, anc4, " | ".join(ancs), ",".join(gids)))
            leaves.append(node)
    print(len(leaves))
    fhi.close()
    fho.close()

if __name__ == "__main__":
    corncyc2tsv()

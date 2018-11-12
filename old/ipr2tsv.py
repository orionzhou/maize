#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description = 'parse interproscan output'
    )
    parser.add_argument(
        'fi', help = 'input file'
    )
    parser.add_argument(
        'fo', help = 'output file (tsv)'
    )
    args = parser.parse_args()
    (fi, fo) = (args.fi, args.fo)

    dic = dict()
    fhi = open(fi, "r")
    for line in fhi:
        ps = line.strip().split("\t")
        aid, md5, plen, aly, sname, sdesc, beg, end, score, status, date = ps[0:11]
        iid = ps[11] if len(ps) >= 12 else ''
        gid = ps[13] if len(ps) >= 14 else ''
        if not aid in dic:
            dic[aid] = [set(), set()]
        iids = iid.split(",")
        for iid in iids:
            if iid:
                dic[aid][0].add(iid)
        gids = gid.split("|")
        for gid in gids:
            if gid:
                dic[aid][1].add(gid)
    fhi.close()

    fho = open(fo, "w")
    for aid, lst in dic.items():
        iids, gids = lst
        fho.write("%s\t%s\n" % (aid, ";".join(gids)))
    fho.close()

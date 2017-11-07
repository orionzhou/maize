#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import argparse
import re

if __name__ == '__main__':
    dirw = '/home/springer/zhoux379/data/misc2/urt'
    fi = op.join(dirw, 'checklater.tsv')
    fp = op.join(dirw, '20.ped.raw.tsv')
    fo = op.join(dirw, 'checklater_out.tsv')
    fhp = open(fp, "r")
    alltxt = fhp.read()
    fhp.close()

    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        ps = line.strip().split("\t")
        sid = ps[0]
        lines = re.findall(r"(^.*?" + re.escape(sid) + ".*?$)", alltxt, re.MULTILINE)
        if len(lines) > 20:
            lines = lines[0:20]
        fho.write("%s" % line)
        fho.write("\n".join(lines) + "\n\n")
    fhi.close()
    fho.close()

 

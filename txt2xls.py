#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import xlwt
import xlrd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert plain txt file to *.xls'
    )
    parser.add_argument('fi', help = 'input file (*.txt)')
    parser.add_argument('fo', help = 'output file (*.xls)')
    args = parser.parse_args()
    fi, fo = args.fi, args.fo
    fhi = open(fi, "r")
    rows = []
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Sheet1')
    i = 0
    for line in fhi:
        line = line.strip("\n")
        ps = line.split("\t")
        for j in range(len(ps)):
            ws.write(i, j, ps[j])
        i += 1
    wb.save(fo)

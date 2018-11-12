#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import openpyxl as xl

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'report fasta sequence lengths'
    )
    parser.add_argument(
        'fi', help = 'input file (*.xlsx)'
    )
    parser.add_argument(
        'fo', help = 'output file (*.xls)'
    )
    args = parser.parse_args()
    fi, fo = args.fi, args.fo
    wb = xl.load_workbook(fi)
    wb.save(fo)

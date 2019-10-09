#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read and write JSON file.
"""

import os.path as op
import sys
import logging
import json

from maize.formats.base import must_open, ndigit, prettysize

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read and write json files'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('fastp', help='Convert fastp output(*.json) files to tsv file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('json', nargs='+', help='one or more (fastp) json file(s)')
    sp1.set_defaults(func = fastp)

    sp1 = sp.add_parser('bbduk', help='Convert bbduk output(*.json) files to tsv file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('json', nargs='+', help='one or more (bbduk) json file(s)')
    sp1.add_argument('--skip', type = int, default = 1, help = 'number of lines to skip')
    sp1.set_defaults(func = bbduk)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

def fastp(args):
    """
    %prog fastp jsonfile

    Convert fastp json to tsv file.
    """
    jsons = args.json
    logging.info("reading %s files..." % len(jsons))
    keys = """passed_filter_reads
        low_quality_reads
        too_many_N_reads
        too_short_reads
        too_long_reads""".split()
    print('\t'.join(['sid'] + keys))
    for fi in jsons:
        sid = op.basename(op.splitext(fi)[0])
        fhi = must_open(fi)
        js = json.load(fhi)
        print("\t".join([sid] + [str(js['filtering_result'][x]) for x in keys]))

def bbduk(args):
    """
    %prog bbduk jsonfile

    Convert bbduk json to tsv file.
    """
    jsons = args.json
    skip = args.skip
    logging.info("reading %s files..." % len(jsons))
    keys = "readsIn readsRemoved readsOut ".split()
    print('\t'.join(['sid'] + keys))
    for fi in jsons:
        sid = op.basename(op.splitext(fi)[0])
        fhi = must_open(fi)
        if skip >= 1:
            for i in range(skip):
                next(fhi)
        js = json.load(fhi)
        print("\t".join([sid] + [str(js[x]) for x in keys]))

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read and write JSON file.
"""

import os.path as op
import sys
import logging
import json
import pandas as pd

from maize.formats.base import must_open, ndigit, prettysize

def main():
    import argparse
    ps = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read and write json files'
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('meta', help='Convert study meta table to IGV track configuration json',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help='input meta table')
    sp1.add_argument('fo', help='output track json')
    sp1.add_argument('--label', default='test', help='name of track sets')
    sp1.add_argument('--desc', default='sample descriptions', help='description of track sets')
    sp1.add_argument('--type', default='annotation', help='track type')
    sp1.add_argument('--format', default='gff3', help='track format')
    sp1.add_argument('--height', default=50, help='track height')
    sp1.add_argument('--displayMode', default='SQUISHED', help='display mode')
    sp1.add_argument('--url_prefix', default="https://s3.msi.umn.edu/zhoup-igv/Zmays-B73", help='S3 URL prefix')
    sp1.add_argument('--url_suffix', default="gff.gz", help='S3 URL suffix')
    sp1.add_argument('--index', action='store_true', help='has index?')
    sp1.add_argument('--idx_suffix', default="tbi", help='index suffix')
    sp1.set_defaults(func = meta2json)

    sp1 = sp.add_parser('fastp', help='Convert fastp output(*.json) files to tsv file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('json', nargs='+', help='one or more (fastp) json file(s)')
    sp1.set_defaults(func = fastp)

    sp1 = sp.add_parser('bbduk', help='Convert bbduk output(*.json) files to tsv file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('json', nargs='+', help='one or more (bbduk) json file(s)')
    sp1.add_argument('--skip', type = int, default = 1, help = 'number of lines to skip')
    sp1.set_defaults(func = bbduk)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()

def meta2json(args):
    """
    %prog meta2json meta json

    Convert meta table to json.
    """
    fi, fo = args.fi, args.fo
    assert op.isfile(fi), "samplelist not found: %s" % fi
    cvts = dict(SampleID=str,Tissue=str,Genotype=str,Treatment=str,Replicate=int,paired=bool,spots=float,avgLength=float)
    sl = pd.read_csv(fi, sep="\t", header=0, converters=cvts)
    od = dict(label=args.label, description=args.desc)
    tracks = []
    for i in range(len(sl)):
        sid, tis, gt, cond, rep = sl['SampleID'][i], sl['Tissue'][i], \
            sl['Genotype'][i], sl['Treatment'][i], sl['Replicate'][i]
        name = "%s %s %s rep%d (%s)" % (tis, gt, cond, rep, sid)
        track = dict(name=name,height=args.height,
                     displayMode=args.displayMode)
        if args.type == 'bigwig2':
            track['type'] = 'merged'
            url1 = "%s/%s.plus.%s" % (args.url_prefix, sid, args.url_suffix)
            url2 = "%s/%s.minus.%s" % (args.url_prefix, sid, args.url_suffix)
            dict1 = dict(type='wig',format='bigwig',url=url1,color="red")
            dict2 = dict(type='wig',format='bigwig',url=url2,color="green")
            track['tracks'] = [dict1, dict2]
        else:
            if args.type == 'bigwig':
                track['type'] = 'wig'
            else:
                track['type'] = args.type
            track['format'] = args.format
            track['url'] = "%s/%s.%s" % (args.url_prefix, sid, args.url_suffix)
            if args.index:
                track['indexURL'] = "%s/%s.%s.%s" % (args.url_prefix, sid, args.url_suffix, args.idx_suffix)
        tracks.append(track)
    od['tracks'] = tracks
    fho = open(fo, "w")
    fho.write(json.dumps(od, indent=2))
    fho.close()

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

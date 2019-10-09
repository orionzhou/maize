#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path as op
import sys
import yaml
import pandas as pd
import primer3

def ufmu(args):
    fhc = open(args.cfg, 'r')
    gcfg = yaml.load(fhc, Loader=yaml.SafeLoader)

    fi = args.fi
    df = pd.read_csv(fi, sep="\t", header=0)
    pkeys = '''
        PRIMER_LEFT_NUM_RETURNED PRIMER_RIGHT_NUM_RETURNED PRIMER_PAIR_NUM_RETURNED
        MER_LEFT_0 PRIMER_RIGHT_0 PRIMER_LEFT_0_SEQUENCE PRIMER_RIGHT_0_SEQUENCE
    '''.split()
    hcols = '''pid start end n_left n_right n_pair
        left.start left.size left.seq right.start right.size right.seq'''.split()
    print("\t".join(hcols))
    for i in range(len(df)):
        pid, start, end, seq = df['pid'][i], int(df['start'][i]), int(df['end'][i]), df['seq'][i]
        scfg = dict()
        scfg['SEQUENCE_ID'] = pid
        scfg['SEQUENCE_TEMPLATE'] = seq
        scfg['SEQUENCE_TARGET'] = [start+1, end-start]
        x = primer3.bindings.designPrimers(scfg, gcfg)
        n_left = x['PRIMER_LEFT_NUM_RETURNED']
        n_right = x['PRIMER_RIGHT_NUM_RETURNED']
        n_pair = x['PRIMER_PAIR_NUM_RETURNED']
        left_start, left_size = x['PRIMER_LEFT_0']
        right_start, right_size = x['PRIMER_RIGHT_0']
        left_seq = x['PRIMER_LEFT_0_SEQUENCE']
        right_seq = x['PRIMER_RIGHT_0_SEQUENCE']
        vs = [pid, start, end,
              n_left, n_right, n_pair,
              left_start-1, left_size, left_seq,
              right_start-1, right_size, right_seq]
        print("\t".join([str(v) for v in vs]))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'primer3 utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("ufmu",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "design primers for UfMu insertions")
    sp1.add_argument('fi', help = 'input sequence table (`pid`, `start.t`, `end.t`, `seq`)')
    sp1.add_argument('--cfg', default='/home/springer/zhoux379/git/maize/apps/primer3.yml', help = 'primer3 global config')
    sp1.set_defaults(func = ufmu)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


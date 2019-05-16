#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
import pandas as pd
from multiprocessing import Pool
import egglib
egglib.wrappers.paths['codeml'] = 'codeml'

def egglib_stat(fi, skeys='S thetaW Pi D lseff nseff'.split()):
    aln = egglib.io.from_fasta(fi, cls=egglib.Align, groups=True)

    cs = egglib.stats.ComputeStats()
    for skey in skeys:
        cs.add_stats(skey)
    stats = cs.process_align(aln)

    try:
        res = egglib.wrappers.codeml(aln, None, 'M0')
        stats['omega'] = res['omega']
    except:
        stats['omega'] = None

    return stats

def stat1(gid, sids, db):
    ftmp = "%s.fas" % gid
    fht = open(ftmp, 'w')
    for sid in sids:
        sgid = '%s#%s' % (sid, gid)
        seq = db[sgid][:-3].seq
        fht.write(">%s\n" % sid)
        fht.write(seq + "\n")
    fht.close()
    stats = egglib_stat(ftmp)
    os.remove(ftmp)
    return stats

def stats(args):
    fi, fg, fs, fo = args.fi, args.fg, args.fs, args.fo
    batch, batch_size = args.batch, args.batch_size
    db = Fasta(fi)

    gids = [line.rstrip('\n') for line in open(fg,'r')]
    sids = [line.rstrip('\n') for line in open(fs,'r')]

    fho = open(fo, 'w')
    skeys = 'lseff nseff S thetaW Pi D omega'.split()
    fho.write("\t".join(['gid'] + skeys) + "\n")

    # pool = Pool(nthreads)
    # inputs = [(gid,sids,db) for gid in gids]
    # for stats in pool.imap_unordered(stat1, inputs):
    i = (batch-1) * batch_size + 1
    j = batch * batch_size
    for idx in range(i-1,j):
        gid = gids[idx]
        stats = stat1(gid, sids, db)
        svals = [str(stats[skey]) for skey in skeys]
        fho.write("\t".join([gid] + svals) + "\n")

    fho.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'popgen utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("stats",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "compute diversity stats for CDS alignment")
    sp1.add_argument('fi', help = 'input fasta db')
    sp1.add_argument('fg', help = 'gene ID list')
    sp1.add_argument('fs', help = 'sample list')
    sp1.add_argument('fo', help = 'output file (*.tsv)')
    sp1.add_argument('--batch', type=int, default=1, help = 'batch to process')
    sp1.add_argument('--batch_size', type=int, default=1000, help = 'number genes per batch')
    sp1.set_defaults(func = stats)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


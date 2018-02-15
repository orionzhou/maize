#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re
from pyfaidx import Fasta

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
def read_bed(fb):
    locs = []
    fhb = open(fb, "r")
    for line in fhb:
        if line.startswith("#"):
            continue
        row = line.strip("\n").split("\t")
        assert len(row) >= 3, "not >=3 fields:\n%s" % line
        rs = row[0:3]
        rs[1] = int(rs[1]) + 1
        rs[2] = int(rs[2])
        if len(row) >= 4 and row[3] != '':
            rs.append(row[3])
        locs.append(rs)
    return locs
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'random retrieve fasta sequences'
    )
    parser.add_argument(
            'db', help = 'sequence database (fasta file or genome ID)'
    )
    parser.add_argument(
            'loc', help = 'location string(s) or BED file(s) (separated by ",")'
    )
    parser.add_argument(
            'out', nargs = "?", help = 'output file [stdout]'
    )
    parser.add_argument(
            '--padding', action = "store_true", help = 'padding to size [No]'
    )
    args = parser.parse_args()

    db = ""
    if op.isfile(args.db):
        db = Fasta(args.db)
    else:
        f_db = "%s/%s/11_genome.fas" % (os.environ["genome"], args.db)
        assert op.isfile(f_db), "cannot find %s" % args.db
        db = Fasta(f_db)

    if args.out and args.out != '-' and args.out != 'stdout':
        sys.stdout = open(args.out, "w")

    reg1 = re.compile("^([\w\-]+)\:([\d,]+)(\-|\.{1,2})([\d,]+)$")
    reg2 = re.compile("^([\w\-]+)$")
    locs = []
    for loc in args.loc.split(","):
        if op.isfile(loc):
            locs += read_bed(loc)
        else:
            res = reg1.match(loc)
            if res:
                sid, beg, end = res.group(1), res.group(2), res.group(4)
                beg = int(beg.replace(",", ""))
                end = int(end.replace(",", ""))
                locs.append([sid, beg, end])
            else:
                res = reg2.match(loc)
                if res:
                    sid = res.group(1)
                    beg = 1
                    if sid in db:
                        end = len(db[sid])
                        locs.append([sid, beg, end])
                    else:
                        eprint("%s not in db => skipped" % sid)
                else:
                    eprint("%s: unknown locstr => skipped" % loc)

    for loc in locs:
        sid, beg, end = loc[0:3]
        oid = "%s-%d-%d" % (sid, beg, end)
        if len(loc) > 3:
            oid = loc[3]
        if sid not in db:
            eprint("%s not in db => skipped" % sid)
            continue
        size = end - beg + 1
        bp_pad = 0
        if beg < 1:
            bp_pad += 1 - beg
            beg = 1
        if beg > len(db[sid]):
            bp_pad = 1
            beg = len(db[sid])
        if end > len(db[sid]):
            bp_pad += end - len(db[sid])
            end = len(db[sid])
        seq = db[sid][beg-1:end].seq
        if args.padding:
            if bp_pad > 0:
                if end-beg+1 < 30:
                    seq = "N" * size
                else:
                    seq += "N" * bp_pad
            assert len(seq) == size, "error in seq size: %s:%d-%d %d" % (sid, beg, end, bp_pad)
        print(">%s" % oid)
        print(seq)

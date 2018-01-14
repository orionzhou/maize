#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'convert inbred vcf to hybrid vcf'
    )
    parser.add_argument(
        'fi', help = 'input file (.vcf)'
    )
    parser.add_argument(
        'fo', help = 'output file (.bed)'
    )
    parser.add_argument(
        '--noindel', action = 'store_true', help = 'remove InDel [No]'
    )
    args = parser.parse_args()

    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith("#"):
            continue
        row = line.strip("\n").split("\t")
        sid, pos, ref, alt = row[0], row[1], row[3], row[4]
        pos = int(pos)
        vnt = "M:%s" % alt
        if len(ref) > 1:
            assert len(alt) == 1, "error: %s" % line
            vnt = "D:%d" % (len(ref) - 1)
            pos = pos + 1
            if args.noindel:
                continue
        elif len(alt) > 1:
            assert len(ref) == 1, "error: %s" % line
            vnt = "I:%s" % alt[1:]
            if args.noindel:
                continue
        else:
            assert len(ref) == 1 and len(alt) == 1, "error: %s" % line
        fho.write("\t".join([sid, str(pos-1), str(pos), vnt]) + "\n")
    fhi.close()
    fho.close()


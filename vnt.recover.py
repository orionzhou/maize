#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse

def add_one_vnt(rdic, pos, ref, alt, qd):
    if pos in rdic:
        pref, palt, pqd = rdic[pos]
        if qd > pqd:
            rdic[pos] = [ref, alt, qd]
    else:
        rdic[pos] = [ref, alt, qd]
def vnt_tile(vntstr, seqid, ibeg):
    rdic = dict()
    if vntstr == '':
        return rdic
    vnts = vntstr.split(" ")
    for vnt in vnts:
        seqid1, beg, end, ref, alt, qd = vnt.split(":")
        if qd == 'NA': qd = '1'
        beg, end, qd = int(beg), int(end), float(qd)
        assert seqid == seqid1, "%s not on %s" % (vntstr, seqid)
        assert end - beg + 1 == len(ref), "%s length not right" % vntstr
        if len(ref) == 1 and len(alt) == 1:
            pos = beg - ibeg + 1
            add_one_vnt(rdic, pos, ref, alt, qd)
        elif len(ref) > 1 and len(alt) == 1:
            for apos in range(beg+1,end+1):
                rpos = apos - ibeg + 1
                add_one_vnt(rdic, rpos, ref[apos-beg], '', qd)
        elif len(ref) == 1 and len(alt) > 1:
            pos = beg - ibeg + 1
            add_one_vnt(rdic, pos, ref, alt, qd)
        else:
            print("unknown mnp: %s" % vntstr)
            sys.exit(1)
    return rdic
def vnt_apply(seq, vdic):
    vseq = list(seq)
    for pos, vnt in vdic.items():
        if pos < 1 or pos > len(seq):
            continue
        ref, alt, qd = vnt
        idx = pos - 1
        assert ref == vseq[idx], "refseq conflict: %d:%s:%s" % (pos, ref, alt)
        vseq[idx] = alt
    return ''.join(vseq)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'recover fasta from ref-seq + variants'
    )
    parser.add_argument(
        'fi', help = 'input file (.tsv)'
    )
    parser.add_argument(
        'fo', help = 'output file (.tsv)'
    )
    args = parser.parse_args()
    fi, fo = args.fi, args.fo

    fhi = open(fi, "r")
    fho = open(fo, "w")
    for line in fhi:
        line = line.strip("\n")
        ary = line.split("\t")
        if len(ary) != 8:
            print("not 8 fields:\n$line")
            sys.exit(1)
        idx, seqid, beg, end, srd, tid, seq, vntstr = ary
        beg, end = int(beg), int(end)
        assert end-beg+1 == len(seq), "seq len error: %s:%s:%s" % (seqid, beg, end)
        vdic = vnt_tile(vntstr, seqid, beg)
        nseq = vnt_apply(seq, vdic)
        fho.write("%s\t%s\n" % (idx, nseq))
    fhi.close()
    fho.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import re
import string

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid
from pyfaidx import Fasta
import pandas as pd

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

FastaExt = ("fasta", "fa", "fas", "fna", "cds", "pep", "faa", "fsa", "seq", "nt", "aa")

def ndigit(num):
    if num < 1:
        print("no digits: %g" % num)
        sys.exit(1)
    digit = 0
    while num >= 1:
        num /= 10.0
        digit += 1
    return digit

def prettysize(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

def size(args):
    if args.header:
        print("seqid\tsize")
    fname, fext = op.splitext(args.fi)
    fho = must_open(args.fo, 'w')
    if args.fi in ['stdin', '-'] or fext in ['.gz','.bz2']:
        fh = must_open(args.fi)
        for rcd in SeqIO.parse(fh, "fasta"):
            sid, size = rcd.id, len(rcd)
            if args.bed:
                fho.write("%s\t%d\t%d\n" % (sid, 0, size))
            else:
                fho.write("%s\t%d\n" % (sid, size))
    elif fext in [".%s" % x for x in FastaExt]:
        from pyfaidx import Fasta
        fas = Fasta(args.fi)
        for sid in fas.keys():
            size = len(fas[sid])
            if args.bed:
                fho.write(f"{sid}\t0\t{size}\n")
            else:
                fho.write(f"{sid}\t{size}\n")
    else:
        logging.error("%s is not a supported format" % fext)

def desc(args):
    fh = must_open(args.fi)
    fho = must_open(args.fo)
    if args.header:
        fho.write("seqid\tdesc\n")
    for rcd in SeqIO.parse(fh, "fasta"):
        sid, desc = rcd.id, rcd.description
        if sid == desc:
            desc = ''
        fho.write("%s\t%s\n" % (sid, desc))

def clean(args):
    import re
    reg = re.compile("[^ATCGN]")
    fhi = must_open(args.fi)
    fho = must_open(args.fo, 'w')
    cnt = 0
    for rcd in SeqIO.parse(fhi, "fasta"):
        sid, seq = rcd.id, str(rcd.seq).upper()
        newseq, ncnt = reg.subn("N", seq)
        cnt += ncnt
        nrcd = SeqRecord(Seq(newseq), id = sid, description = "")
        SeqIO.write(nrcd, fho, "fasta")
    logging.debug("Total bad char: %d" % cnt)

def translate(args):
    fh = must_open(args.fi)
    for rcd in SeqIO.parse(fh, "fasta"):
        sid = rcd.id
        aa = rcd.seq.translate(to_stop = True)
        nrcd = SeqRecord(aa, id = sid, description = "")
        SeqIO.write(nrcd, sys.stdout, "fasta")

def extract(args):
    import re
    from jcvi.formats.bed import Bed

    db = ""
    if op.isfile(args.db):
        db = Fasta(args.db)
    else:
        f_db = "%s/data/%s/10_genome.fna" % (os.environ["genome"], args.db)
        assert op.isfile(f_db), "cannot find %s" % args.db
        db = Fasta(f_db)

    reg1 = re.compile("^([\w\-]+)\:([\d,]+)(\-|\.{1,2})([\d,]+)$")
    reg2 = re.compile("^([\w\-]+)$")
    bed = Bed()
    if op.isfile(args.loc):
        if args.list:
            fho = must_open(args.loc, 'r')
            for line in fho:
                sid = line.strip()
                beg = 0
                if sid in db:
                    end = len(db[sid])
                    bed.add("%s\t%d\t%d\n" % (sid, beg, end))
                else:
                    logging.warning(f"{sid} not in db => skipped")
        else:
            bed = Bed(args.loc, sorted=False)
    else:
        for loc in args.loc.split(","):
            res = reg1.match(loc)
            if res:
                sid, beg, end = res.group(1), res.group(2), res.group(4)
                beg = int(beg.replace(",", ""))
                end = int(end.replace(",", ""))
                bed.add("%s\t%d\t%d\n" % (sid, beg-1, end))
            else:
                res = reg2.match(loc)
                if res:
                    sid = res.group(1)
                    beg = 0
                    if sid in db:
                        end = len(db[sid])
                        bed.add("%s\t%d\t%d\n" % (sid, beg, end))
                    else:
                        logging.warning(f"{sid} not in db => skipped")
                else:
                    logging.warning(f"{loc}: unknown locstr => skipped")

    rcds = []
    for b in bed:
        sid, beg, end = b.seqid, b.start, b.end
        oid = sid if args.list else f"{sid}-{beg}-{end}"
        if b.accn:
            oid = b.accn
        if sid not in db:
            logging.warning(f"{sid} not in db => skipped")
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

        if args.tsv:
            print("\t".join([sid, str(beg), str(end), seq]))
        else:
            rcd = SeqRecord(Seq(seq), id = oid, description = '')
            SeqIO.write([rcd], sys.stdout, 'fasta')

def split_old(args):
    fi, dirw = op.realpath(args.fi), op.realpath(args.outdir)
    n = args.N
    if not op.exists(dirw):
        mkdir(dirw)
    else:
        sh("rm -rf %s/*" % dirw)

    cdir = os.path.dirname(os.path.realpath(__file__))
    cwd = os.getcwd()
    os.chdir(dirw)

    sh("ln -sf %s part.fas" % fi)
    sh("pyfasta split -n %d part.fas" % n)
    sh("rm part.fas part.fas.*")

    digit = ndigit(n)
    sizes = []
    for i in range(0,n):
        fmt = "part.%%0%dd.fas" % digit
        fp = fmt % i
        sizes.append(os.stat(fp).st_size)
    sizes.sort()
    print("size range: %s - %s" % (prettysize(sizes[0]), prettysize(sizes[n-1])))

def tile(args):
    from jcvi.utils.location import maketile

    fhi = must_open(args.fi)
    winstep, winsize = args.step, args.size
    for rcd in SeqIO.parse(fhi, "fasta") :
        size = len(rcd.seq)
        sid, beg, end = rcd.id, 1, size
        ary = rcd.id.split("-")
        if len(ary) >= 3:
            sid, beg, end = ary[0], int(ary[1]), int(ary[2])
            assert size == end - beg + 1, "size error: %s not %d" % (rcd.id, size)
        elif len(ary) == 2:
            sid, beg = ary[0], int(ary[1])
            end = beg + size - 1

        wins = maketile(1, size, winsize, winstep)
        rcds = []
        seq = str(rcd.seq)
        for rbeg, rend in wins:
            abeg, aend = beg + rbeg - 1, beg + rend - 1
            ssid = "%s-%d-%d" % (sid, abeg, aend)
            seqstr = seq[rbeg-1:rend]
            rcds.append(SeqRecord(Seq(seqstr), id = ssid, description = ''))
        SeqIO.write(rcds, sys.stdout, "fasta")
    fhi.close()

def merge(args):
    cfg = args.cfg
    for line in must_open(cfg):
        line = line.strip(" \t\n\r")
        if line == "":
            continue
        (pre, fseq) = line.split(",")
        if not os.access(fseq, os.R_OK):
            print("no access to input file: %s" % fseq)
            sys.exit(1)

        fh = must_open(fseq)
        seq_it = SeqIO.parse(fh, "fasta")
        seqs = [SeqRecord(rcd.seq, id = pre + "|" + rcd.id,
            description = '') for rcd in seq_it]
        SeqIO.write(seqs, sys.stdout, "fasta")
        fh.close()

def merge_pe(args):
    (fi1, fi2, fo) = (args.fi1, args.fi2, args.fo)
    assert op.isfile(fi1), "cannot read %s" % fi1
    assert op.isfile(fi2), "cannot read %s" % fi2

    fhi2 = open(fi1, "rb")
    fhi1 = open(fi2, "rb")
    fho = open(fo, "wb")

    for lst1, lst2 in zip(read_fasta(fi1, fhi1), read_fasta(fi2, fhi2)):
        seqid1, seq1 = lst1
        seqid2, seq2 = lst2
        assert seqid1 == seqid2 and len(seq1) == len(seq2), \
                "%s: seq[%d] != %s: seq[%d]" % \
                (seqid1, len(seq1), seqid2, len(seq2))
        if args.join:
            fho.write((">%s\n%s\n" % (seqid1, seq1+seq2)).encode('utf8'))
        else:
            if args.nosuffix:
                fho.write((">%s\n%s\n" % (seqid1, seq1)).encode('utf8'))
                fho.write((">%s\n%s\n" % (seqid2, seq2)).encode('utf8'))
            else:
                fho.write((">%s.1\n%s\n" % (seqid1, seq1)).encode('utf8'))
                fho.write((">%s.2\n%s\n" % (seqid2, seq2)).encode('utf8'))
    fhi1.close()
    fhi2.close()
    fho.close()

def gaps(args):
    import re
    reg = re.compile("N+")
    fh = must_open(args.fi)
    fho = must_open(args.fo, 'w')
    for rcd in SeqIO.parse(fh, "fasta"):
        sid, seq = rcd.id, str(rcd.seq).upper()
        for res in reg.finditer(seq):
            beg, end = res.start(0), res.end(0)
            if end - beg >= args.gap:
                fho.write("%s\t%d\t%d\n" % (sid, beg, end))

def fas2aln(args):
    from Bio import AlignIO
    fh = must_open(args.fi)
    alns = AlignIO.parse(fh, "fasta")
    AlignIO.write(alns, sys.stdout, "clustal")
    fh.close()

def extract_chrom_num(sid, opt):
    chrom = False
    if opt in 'Brapa Ppersica'.split():
        ptn = "^[AG]([0-9]{1,2})"
        res = re.search(ptn, sid, re.IGNORECASE)
        chrom = res.group(1).lstrip('0') if res else False
    elif opt == 'Sitalica':
        ptn = "^([IVX]{1,4})"
        dic_chrom = dict(I=1,II=2,III=3,IV=4,V=5,VI=6,VII=7,VIII=8,IX=9)
        res = re.search(ptn, sid, re.IGNORECASE)
        chrom = str(dic_chrom[res.group(1)]) if res else False
    elif opt == 'Vvinifera':
        ptn = "^([1-9][0-9]{0,1})$"
        res = re.search(ptn, sid, re.IGNORECASE)
        chrom = res.group(1) if res else False
    else:
        ptn = "^(chr|chromsome|Ta)?[ _]*(0*[1-9XY][0-9]{0,1}[A-Z]?)" #(MtrunA17)?
        res = re.search(ptn, sid, re.IGNORECASE)
        chrom = res.group(2).lstrip('0') if res else False
    return chrom

def rename(args):
    import re
    from pyfaidx import Fasta

    fi, fo, fmf, fmb = args.fi, args.fo, args.fmf, args.fmb
    opt, merge_short, gap = args.opt, args.merge_short, args.gap
    prefix_chr, prefix_ctg = args.prefix_chr, args.prefix_ctg

    db = Fasta(fi)

    sdic, cdic = dict(), dict()
    ccnt = 1
    for sid in db.keys():
        size = len(db[sid])
        chrom = extract_chrom_num(sid, opt)
        if chrom:
            num = chrom.strip(string.ascii_letters)
            num = int(num) if len(num) > 0 else 0
            sdic[sid] = [chrom, size, num]
        else:
            cdic[sid] = [ccnt, size]
            ccnt += 1

    maxnum = 0
    if len(sdic) > 0:
        maxnum = max(v[2] for k,v in sdic.items())
        assert maxnum <= 99, ">99 [%d] chroms: not supported" % maxnum
        for sid, sval in sdic.items():
            chrom, size, num = sval
            sdic[sid][0] = f"{prefix_chr}0{chrom}" if num <= 9 and maxnum >= 10 else f"{prefix_chr}{chrom}"

    slst = sorted(sdic.items(), key = lambda t: t[1][0])
    clst = sorted(cdic.items(), key = lambda t: t[1][0])
    cdigits = ndigit(clst[-1][1][0]) if len(clst) > 0 else 1
    cfmt = "%s%%0%dd" % (prefix_ctg, cdigits)
    logging.debug("%d chromosomes, %d scaffolds/contigs" % (len(sdic), len(cdic)))

    fname, fext = op.splitext(fi)
    if fext not in [".%s" % x for x in FastaExt]:
        logging.error("%s is not a supported format" % fext)
        sys.exit(1)

    fho = open(fo, "w")
    fhf = open(fmf, "w")
    fhb = open(fmb, "w")

    i = 1
    if len(slst) > 0:
        for sid, sval in slst:
            nsid, size, num = sval
            fhf.write(f"{sid}\t0\t{size}\t+\t{nsid}\t0\t{size}\t{i}\n")
            fhb.write(f"{nsid}\t0\t{size}\t+\t{sid}\t0\t{size}\t{i}\n")
            nrcd = SeqRecord(Seq(str(db[sid])), id = nsid, description = '')
            SeqIO.write(nrcd, fho, "fasta")
            i += 1

    if len(clst) > 0 and merge_short:
        zid = f"{prefix_chr}x" if maxnum <= 9 else f"{prefix_chr}99"
        pos = 0
        seq = ''
        for cid, sval in clst:
            ccnt, size = sval
            start, end = pos, pos + size
            if pos > 0:
                start += gap
                end += gap
                seq += "N" * gap
            seq += str(db[cid])
            fhf.write("%s\t%d\t%d\t+\t%s\t%d\t%d\t%d\n" % (cid, 0, size, zid, start, end, i))
            fhb.write("%s\t%d\t%d\t+\t%s\t%d\t%d\t%d\n" % (zid, start, end, cid, 0, size, i))
            pos = end
            i += 1
        nrcd = SeqRecord(Seq(seq), id = zid, description = '')
        SeqIO.write(nrcd, fho, "fasta")
    else:
        for cid, sval in clst:
            ccnt, size = sval
            ncid = cfmt % ccnt
            fhf.write("%s\t%d\t%d\t+\t%s\t%d\t%d\t%d\n" % (cid, 0, size, ncid, 0, size, i))
            fhb.write("%s\t%d\t%d\t+\t%s\t%d\t%d\t%d\n" % (ncid, 0, size, cid, 0, size, i))
            nrcd = SeqRecord(Seq(str(db[cid])), id = ncid, description = '')
            SeqIO.write(nrcd, fho, "fasta")
            i += 1

    fhf.close()
    fhb.close()
    fho.close()

def rmdot(args):
    from string import maketrans
    fh = must_open(args.fi)
    tt = maketrans(".", "-")
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            print(line)
        else:
            print(line.translate(tt))
    fh.close()

def cleanid(args):
    fh = must_open(args.fi)
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            print(re.sub('[(:.].*$', '', line))
        else:
            print(line)
    fh.close()

def suffix(args):
    fh = must_open(args.fi)
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            print(line + args.suf)
        else:
            print(line)
    fh.close()

def dedup(args):
    fh = must_open(args.fi)
    seqd = dict()
    seq_it = SeqIO.parse(fh, "fasta")
    cnt, cntu = 0, 0
    for rcd in seq_it:
        rid, seq = rcd.id, rcd.seq
        cnt += 1
        if seq in seqd:
            seqd[seq] += args.sep + rid
        else:
            seqd[seq] = rid
            cntu += 1
    fh.close()

    seqs = [SeqRecord(seq, id=rids, description='') for seq, rids in seqd.items()]
    SeqIO.write(seqs, sys.stdout, "fasta")
    logging.info("%6d total sequences" % cnt)
    logging.info("%6d non-redundant sequences" % cntu)

def add_nr(args):
    fh = must_open(args.fi)
    seqd = dict()
    seq_it = SeqIO.parse(fh, "fasta")
    cnt, cntu = 0, 0
    for rcd in seq_it:
        rid, seq = rcd.id, rcd.seq
        cnt += 1
        if seq in seqd:
            seqd[seq] += args.sep + rid
        else:
            seqd[seq] = rid
            cntu += 1
    fh.close()

    cntn = 0
    fh = must_open(args.fi2)
    seq_it = SeqIO.parse(fh, "fasta")
    for rcd in seq_it:
        rid, seq = rcd.id, rcd.seq
        if seq not in seqd:
            seqd[seq] = rid
            cntn += 1
    fh.close()

    seqs = [SeqRecord(seq, id=rids, description='') for seq, rids in seqd.items()]
    SeqIO.write(seqs, sys.stdout, "fasta")
    logging.info("%6d total sequences" % cnt)
    logging.info("%6d non-redundant sequences" % cntu)
    logging.info("%6d new sequences added" % cntn)

def cds2gene(args):
    fi, fg, fo = args.fi, args.fg, args.fo
    pre = args.prefix
    db = Fasta(fi)
    tg = pd.read_csv(fg, sep='\t', header=0)

    sdic = dict()
    fho = open(fo, "w")
    for i in range(len(tg)):
        if tg['etype'][i] != 'CDS': continue
        gid, chrom, start, end, srd = [tg[col][i] for col in 'gid chrom start end srd'.split()]
        loc = "%s:%d-%d" % (chrom, start, end)
        if gid not in sdic: sdic[gid] = {}
        if 'srd' not in sdic[gid]: sdic[gid]['srd'] = srd
        if 'loc' not in sdic[gid]: sdic[gid]['loc'] = {}
        sdic[gid]['loc'][start] = loc

    print("writing CDS sequences for %d genes" % len(list(sdic.keys())))
    for gid, gdic in sdic.items():
        srd, ldic = gdic['srd'], gdic['loc']
        seqs = []
        if srd == '-':
            locs = sorted(ldic.items(), key=lambda kv: kv[1], reverse=True)
            seqs = [db[loc[1]][:].reverse.complement.seq for loc in locs]
        else:
            locs = sorted(ldic.items(), key=lambda kv: kv[1])
            seqs = [db[loc[1]][:].seq for loc in locs]
        seq = ''.join(seqs)
        rcd_id = gid
        if pre != None and pre != '':
            rcd_id = '%s%s' % (pre, gid)
        nrcd = SeqRecord(Seq(seq), id = rcd_id, description = "")
        SeqIO.write(nrcd, fho, "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fasta utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("size", help = "Report length for each sequence")
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.add_argument('fo', help = 'output file (tab)')
    sp1.add_argument('--header', action = 'store_true', help = 'add header')
    sp1.add_argument('--bed', action = 'store_true', help = 'output in BED')
    sp1.set_defaults(func = size)

    sp1 = sp.add_parser("desc", help = "Report description for each sequence")
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.add_argument('fo', help = 'output file (tsv)')
    sp1.add_argument('--header', action = 'store_true', help = 'add header')
    sp1.set_defaults(func = desc)

    sp1 = sp.add_parser("clean", help = "Remove irregular chararacters")
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.add_argument('fo', help = 'output file (fasta)')
    sp1.set_defaults(func = clean)

    sp1 = sp.add_parser("extract", help = 'retrieve fasta sequences',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('db', help = 'sequence database (fasta or genome ID)')
    sp1.add_argument('loc', help = 'location string(s) or BED file(s) (separated by ",")')
    sp1.add_argument('--padding', action = "store_true", help = 'padding to size')
    sp1.add_argument('--tsv', action = "store_true", help = 'output in tabular format')
    sp1.add_argument('--list', action = "store_true", help = 'input is text file with sequence IDs')
    sp1.set_defaults(func = extract)

    sp1 = sp.add_parser("tile", help = 'create sliding windows that tile the entire sequence',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input fasta file')
    sp1.add_argument('--size', type = int, default = 100, help = 'window size')
    sp1.add_argument('--step', type = int, default = 50, help = 'window step')
    sp1.set_defaults(func = tile)

    sp1 = sp.add_parser("gaps", help = "report gap ('N's) locations in fasta sequences",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.add_argument('fo', help = 'output file (bed)')
    sp1.add_argument('--gap', type = int, default = 10, help = 'min gap size')
    sp1.set_defaults(func = gaps)

    sp1 = sp.add_parser("rename", help = 'rename/normalize sequence IDs, merge short scaffolds/contigs',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input fasta file')
    sp1.add_argument('fo', help = 'output (renamed) fasta file')
    sp1.add_argument('fmf', help = 'forward (old -> new) mapping table of sequence IDs')
    sp1.add_argument('fmb', help = 'backward (new -> old) mapping table of sequence IDs')
    sp1.add_argument('--opt', default = 'default', help = 'chrom renaming optiton')
    sp1.add_argument('--merge_short', action = 'store_true', help = 'merge short scaffolds/contigs')
    sp1.add_argument('--gap', type = int, default = 10000, help = 'number of \'N\'s between short scaffolds/contigs')
    sp1.add_argument('--prefix_chr', default = 'chr', help = 'prefix for renamed sequence IDs')
    sp1.add_argument('--prefix_ctg', default = 'scf', help = 'prefix for short scaffolds/contigs')
    sp1.set_defaults(func = rename)

    sp1 = sp.add_parser("merge", help='merge multiple fasta files and update IDs')
    sp1.add_argument('cfg', help='config file (a text file with identifier followed by the absolute path of fasta in each line)')
    sp1.set_defaults(func = merge)

    sp1 = sp.add_parser("merge_pe",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'merge two (paired) fasta files into one')
    sp1.add_argument('fi1', help='input fasta 1')
    sp1.add_argument('fi2', help='input fasta 2')
    sp1.add_argument('fo', help='output fasta file')
    sp1.add_argument('--nosuffix', action="store_true", help='disable adding *.1 and *.2 sufix')
    sp1.add_argument('--join', action="store_true", help='join seqs of two reads to make one long read')
    sp1.set_defaults(func = merge_pe)

    sp1 = sp.add_parser("rmdot", help='replace periods (.) in an alignment fasta by dashes (-)')
    sp1.add_argument('fi', help='input fasta file')
    sp1.set_defaults(func = rmdot)

    sp1 = sp.add_parser("suffix", help='add a suffix to all sequence identifiers',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help='input fasta file')
    sp1.add_argument('--suf', default='m', help='suffix')
    sp1.set_defaults(func = suffix)

    sp2 = sp.add_parser("cleanid", help='clean sequence IDs in a fasta file')
    sp2.add_argument('fi', help='input fasta file')
    sp2.set_defaults(func = cleanid)

    sp1 = sp.add_parser("2aln", help='convert fasta alignment file to clustal format')
    sp1.add_argument('fi', help='input alignment (.fas)')
    sp1.set_defaults(func = fas2aln)

    sp1 = sp.add_parser("translate", help='translate nucleotide seqs to amino acid seqs')
    sp1.add_argument('fi', help='input fasta file')
    sp1.set_defaults(func = translate)

    sp1 = sp.add_parser("dedup", help='collapse duplicate/identical sequences to a single entry',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help='input fasta file')
    sp1.add_argument('--sep', default=',', help='separator')
    sp1.set_defaults(func = dedup)

    sp1 = sp.add_parser("add_nr", help='add non-redundant sequences from file2 to file1',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help='input fasta file')
    sp1.add_argument('fi2', help='input fasta file 2')
    sp1.add_argument('--sep', default=',', help='separator')
    sp1.set_defaults(func = add_nr)

    sp1 = sp.add_parser("cds2gene", help='concatenate CDS segmetns for each gene')
    sp1.add_argument('fi', help='input fasta file')
    sp1.add_argument('fg', help='gene interval file (*.tsv)')
    sp1.add_argument('fo', help='output fasta file')
    sp1.add_argument('--prefix', default=None, help='prefix added to each record')
    sp1.set_defaults(func = cds2gene)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


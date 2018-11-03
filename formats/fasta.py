#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open, ndigit, prettysize

FastaExt = ("fasta", "fa", "fas", "fna", "cds", "pep", "faa", "fsa", "seq", "nt", "aa")

def size(args):
    if args.header:
        print("seqid\tsize")
    fname, fext = op.splitext(args.fi)
    if args.fi in ['stdin', '-'] or fext in ['.gz','.bz2']:
        fh = must_open(args.fi)
        for rcd in SeqIO.parse(fh, "fasta"):
            sid, size = rcd.id, len(rcd)
            if args.bed:
                print("%s\t%d\t%d" % (sid, 0, size))
            else:
                print("%s\t%d" % (sid, size))
    elif fext in [".%s" % x for x in FastaExt]:
        from pyfaidx import Fasta
        fas = Fasta(args.fi)
        for sid in fas.keys():
            size = len(fas[sid])
            if args.bed:
                print("%s\t%d\t%d" % (sid, 0, size))
            else:
                print("%s\t%d" % (sid, size))
    else:
        logging.error("%s is not a supported format" % fext)

def desc(args):
    fh = must_open(args.fi)
    if args.header:
        print("seqid\tdesc")
    for rcd in SeqIO.parse(fh, "fasta"):
        sid, desc = rcd.id, rcd.description
        if sid == desc:
            desc = ''
        print("%s\t%s" % (sid, desc))

def clean(args):
    import re
    reg = re.compile("[^ATCGN]")
    fh = must_open(args.fi)
    cnt = 0
    for rcd in SeqIO.parse(fh, "fasta"):
        sid, seq = rcd.id, str(rcd.seq).upper()
        newseq, ncnt = reg.subn("N", seq)
        cnt += ncnt
        nrcd = SeqRecord(Seq(newseq), id = sid, description = "")
        SeqIO.write(nrcd, sys.stdout, "fasta")
    logging.debug("Total bad char: %d" % cnt)

def extract(args):
    from pyfaidx import Fasta
    import re
    from maize.formats.bed import Bed
    
    db = ""
    if op.isfile(args.db):
        db = Fasta(args.db)
    else:
        f_db = "%s/%s/10_genome.fna" % (os.environ["genome"], args.db)
        assert op.isfile(f_db), "cannot find %s" % args.db
        db = Fasta(f_db)

    reg1 = re.compile("^([\w\-]+)\:([\d,]+)(\-|\.{1,2})([\d,]+)$")
    reg2 = re.compile("^([\w\-]+)$")
    bed = Bed()
    if op.isfile(args.loc):
        bed = Bed(args.loc)
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
                    beg = 1
                    if sid in db:
                        end = len(db[sid])
                        bed.add("%s\t%d\t%d\n" % (sid, beg, end))
                    else:
                        eprint("%s not in db => skipped" % sid)
                else:
                    eprint("%s: unknown locstr => skipped" % loc)
    
    rcds = []
    for b in bed:
        sid, beg, end = b.seqid, b.start, b.end
        oid = "%s-%d-%d" % (sid, beg, end)
        if b.accn:
            oid = b.accn
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
        rcd = SeqRecord(Seq(seq), id = oid, description = '')
        rcds.append(rcd)
    SeqIO.write(rcds, sys.stdout, 'fasta')

def split(args):
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
    from maize.utils.location import maketile

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
            eprint("no access to input file: %s" % fseq)
            sys.exit(1)
        
        fh = must_open(fseq)
        seq_it = SeqIO.parse(fh, "fasta")
        seqs = [SeqRecord(rcd.seq, id = pre + "|" + rcd.id,
            description = '') for rcd in seq_it]
        SeqIO.write(seqs, sys.stdout, "fasta")
        fh.close()

def gaps(args):
    import re
    reg = re.compile("N+")
    fh = must_open(args.fi)
    for rcd in SeqIO.parse(fh, "fasta"):
        sid, seq = rcd.id, str(rcd.seq).upper()
        for res in reg.finditer(seq):
            beg, end = res.start(0), res.end(0)
            if end - beg >= args.gap:
                print("%s\t%d\t%d" % (sid, beg, end))

def fas2aln(args):
    from Bio import AlignIO
    fh = must_open(args.fi)
    alns = AlignIO.parse(fh, "fasta")
    AlignIO.write(alns, sys.stdout, "clustal")
    fh.close()

def rename(args):
    import re
    from pyfaidx import Fasta
    fi, fs, fo, fmf, fmb = args.fi, args.fs, args.fo, args.fmf, args.fmb
    merge_short, gap = args.merge_short, args.gap
    prefix_chr, prefix_ctg = args.prefix_chr, args.prefix_ctg
    reg_chr = re.compile("^(chr)?([0-9]{1,2})$")
    reg_gb = re.compile("^CM[0-9]+\.[0-9]{1,3}$")
    
    sdic, cdic = dict(), dict()
    scnt, ccnt = 1, 1
    fhs = open(fs, "r")
    for line in fhs:
        sid, size = line.strip().split("\t")
        size = int(size)
        res = reg_chr.match(sid)
        res2 = reg_gb.match(sid)
        if res:
            sdic[sid] = [int(res.group(2)), size]
        elif res2:
            sdic[sid] = [scnt, size]
            scnt += 1
        else:
            cdic[sid] = [ccnt, size]
            ccnt += 1
    slst = sorted(sdic.items(), key = lambda t: t[1][0])
    clst = sorted(cdic.items(), key = lambda t: t[1][0])

    nchrom = slst[-1][1][0]
    sdigits = ndigit(slst[-1][1][0])
    cdigits = ndigit(clst[-1][1][0])
    sfmt = "%s%%0%dd" % (prefix_chr, sdigits)
    cfmt = "%s%%0%dd" % (prefix_ctg, cdigits)
    logging.debug("%d chromosomes, %d scaffolds/contigs" % (len(sdic), len(cdic)))

    fname, fext = op.splitext(fi)
    if fext not in [".%s" % x for x in FastaExt]:
        logging.error("%s is not a supported format" % fext)
        sys.exit(1)
    fas = Fasta(fi)
 
    fho = open(fo, "w")
    fhf = open(fmf, "w")
    fhb = open(fmb, "w")
    for sid, sval in slst:
        scnt, size = sval
        nsid = sfmt % scnt
        fhf.write("%s\t%d\t%d\t+\t%s\t%d\t%d\t%d\n" % (sid, 0, size, nsid, 0, size, scnt))
        fhb.write("%s\t%d\t%d\t+\t%s\t%d\t%d\t%d\n" % (nsid, 0, size, sid, 0, size, scnt))
        nrcd = SeqRecord(Seq(str(fas[sid])), id = nsid, description = '')
        SeqIO.write(nrcd, fho, "fasta")
    i = nchrom + 1
    if merge_short:
        zid = "%sx" % prefix_chr
        if sdigits == 2:
            zid = "%s99" % prefix_chr
        else:
            assert sdigits == 1, "wrong number of chroms: %d" % sdigits
        pos = 0
        seq = ''
        for cid, sval in clst:
            ccnt, size = sval
            start, end = pos, pos + size
            if pos > 0:
                start += gap
                end += gap 
                seq += "N" * gap
            seq += str(fas[cid])
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
            nrcd = SeqRecord(Seq(str(fas[cid])), id = ncid, description = '')
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
            print(line.rstrip(":."))
        else:
            print(line)
    fh.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fasta utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("size", help = "Report length for each sequence") 
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.add_argument('--header', action = 'store_true', help = 'add header')
    sp1.add_argument('--bed', action = 'store_true', help = 'output in BED')
    sp1.set_defaults(func = size)

    sp1 = sp.add_parser("desc", help = "Report description for each sequence") 
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.add_argument('--header', action = 'store_true', help = 'add header')
    sp1.set_defaults(func = desc)

    sp1 = sp.add_parser("clean", help = "Remove irregular chararacters") 
    sp1.add_argument('fi', help = 'input file (fasta)')
    sp1.set_defaults(func = clean)
    
    sp2 = sp.add_parser("extract", 
            help = 'retrieve fasta sequences',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp2.add_argument('db', help = 'sequence database (fasta or genome ID)')
    sp2.add_argument('loc', help = 'location string(s) or BED file(s) (separated by ",")')
    sp2.add_argument('--padding', action = "store_true", help = 'padding to size')
    sp2.set_defaults(func = extract)

    sp3 = sp.add_parser("split", 
            help = 'run pyfasta to split a set of fasta records evenly',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp3.add_argument('fi', help = 'input file (fasta)')
    sp3.add_argument('outdir', help = 'output directory')
    sp3.add_argument('--N', type = int, default = 24, help = 'number pieces')
    #nproc = int(os.environ['nproc'])
    sp3.set_defaults(func = split)

    sp3 = sp.add_parser("tile", 
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'create sliding windows that tile the entire sequence'
    )
    sp3.add_argument('fi', help = 'input fasta file')
    sp3.add_argument('--size', type = int, default = 100, help = 'window size')
    sp3.add_argument('--step', type = int, default = 50, help = 'window step')
    sp3.set_defaults(func = tile)

    sp3 = sp.add_parser("merge", help = 'merge multiple fasta files and update IDs')
    sp3.add_argument('cfg', help = 'config file (a text file with identifier followed by the absolute path of fasta in each line)')
    sp3.set_defaults(func = merge)
 
    sp9 = sp.add_parser("gaps",
            help = "report gap ('N's) locations in fasta sequences",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp9.add_argument('fi', help = 'input file (fasta)')
    sp9.add_argument('--gap', type = int, default = 10, help = 'min gap size')
    sp9.set_defaults(func = gaps)

    sp31 = sp.add_parser("rename",
            help = 'rename/normalize sequence IDs, merge short scaffolds/contigs',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp31.add_argument('fi', help = 'input fasta file')
    sp31.add_argument('fs', help = 'fasta size file (must be generated by "fasta size")')
    sp31.add_argument('fo', help = 'output (renamed) fasta file')
    sp31.add_argument('fmf', help = 'forward (old -> new) mapping table of sequence IDs')
    sp31.add_argument('fmb', help = 'backward (new -> old) mapping table of sequence IDs')
    sp31.add_argument('--merge_short', action = 'store_true', help = 'merge short scaffolds/contigs')
    sp31.add_argument('--gap', type = int, default = 10000, help = 'number of \'N\'s between short scaffolds/contigs')
    sp31.add_argument('--prefix_chr', default = 'chr', help = 'prefix for renamed sequence IDs')
    sp31.add_argument('--prefix_ctg', default = 'ctg', help = 'prefix for short scaffolds/contigs')
    sp31.set_defaults(func = rename)
    
    sp31 = sp.add_parser("rmdot", help = 'replace periods (.) in an alignment fasta by dashes (-)')
    sp31.add_argument('fi', help = 'input fasta file')
    sp31.set_defaults(func = rmdot)
    
    sp32 = sp.add_parser("cleanid", help = 'clean sequence IDs in a fasta file')
    sp32.add_argument('fi', help = 'input fasta file')
    sp32.set_defaults(func = cleanid)
    
    sp11 = sp.add_parser("fas2aln", help = 'convert fasta alignment file to clustal format')
    sp11.add_argument('fi', help = 'input alignment (.fas)')
    sp11.set_defaults(func = fas2aln)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


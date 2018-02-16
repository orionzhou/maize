#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open, digitof_number, sizeof_fmt

def size(args):
    from pyfaidx import Fasta
    fas = Fasta(args.fi)
    if args.header:
        print("seqid\tsize")
    for sid in fas.keys():
        size = len(fas[sid])
        if args.bed:
            print("%s\t%d\t%d" % (sid, 0, size))
        else:
            print("%s\t%d" % (sid, size))

def desc(args):
    fh = must_open(args.fi)
    if args.header:
        print("seqid\tdesc")
    for seqrcd in SeqIO.parse(fh, "fasta"):
        sid, desc = seqrcd.id, seqrcd.description
        if sid == desc:
            desc = ''
        print("%s\t%s" % (sid, desc))

def extract(args):
    from pyfaidx import Fasta
    import re
    from maize.formats.bed import Bed
    
    db = ""
    if op.isfile(args.db):
        db = Fasta(args.db)
    else:
        f_db = "%s/%s/11_genome.fas" % (os.environ["genome"], args.db)
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
        print(">%s" % oid)
        print(seq)

def split(args):
    fi, dirw = op.realpath(args.fi), op.realpath(args.outdir)
    n = args.N
    if not op.exists(dirw):
        makedir(dirw)
    else:
        sh("rm -rf %s/*" % dirw)
    
    cdir = os.path.dirname(os.path.realpath(__file__))
    cwd = os.getcwd()
    os.chdir(dirw)

    sh("ln -sf %s part.fas" % fi)
    sh("pyfasta split -n %d part.fas" % n)
    sh("rm part.fas part.fas.*")

    digit = digitof_number(n)
    sizes = []
    for i in range(0,n):
        fmt = "part.%%0%dd.fas" % digit
        fp = fmt % i
        sizes.append(os.stat(fp).st_size)
    sizes.sort()
    print("size range: %s - %s" % (sizeof_fmt(sizes[0]), sizeof_fmt(sizes[n-1])))

def tile(args):
    fhi = must_open(args.fi)
    winstep, winsize = args.winstep, args.winsize
    for seq in SeqIO.parse(fhi, "fasta") :
        size = len(seq.seq)
        if(float(size) / winstep > 1.3) :
            ary = seq.id.split("-")
            [id, bbeg] = [ary[0], int(ary[1])]

            seqstr = str(seq.seq)
            nf = int(math.ceil(float(size) / piecesize))
            rcds = []
            for i in range(0, nf) :
                rbeg = i * piecesize
                rend = min((i+1) * piecesize, size)
                sseqstr = seqstr[rbeg:rend]
                sid = "%s-%d-%d" % (id, bbeg+rbeg, bbeg+rend-1)
                rcd = SeqRecord(Seq(sseqstr), id = sid, description = '')
                rcds.append(rcd)
                #print "      " + sid
            SeqIO.write(rcds, sys.stdout, "fasta")
        else:
            SeqIO.write(seq, sys.stdout, "fasta")
    fhi.close()

   

def splitlong(args):
    ary = []
    fhi = must_open(args.fi)
    for seq in SeqIO.parse(fhi, "fasta") :
        ary.append(len(seq.seq))
    fhi.close()
    totalsize = sum(ary)
    
    if args.mode == 1:
        piecesize = 100000
    else:
        npieces = 10 * 24 
        piecesize = int((totalsize / npieces) / 100000) * 100000
    print("  total size: %d, size per piece: %d" % (totalsize, piecesize))
    
    fhi = must_open(args.fi)
    fho = open(args.fo, "w")
    for seq in SeqIO.parse(fhi, "fasta") :
        size = len(seq.seq)
        if(float(size) / piecesize > 1.3) :
            print("    splitting %s: %d" %(seq.id, size))
            ary = seq.id.split("-")
            [id, bbeg] = [ary[0], int(ary[1])]

            seqstr = str(seq.seq)
            nf = int(math.ceil(float(size) / piecesize))
            rcds = []
            for i in range(0, nf) :
                rbeg = i * piecesize
                rend = min((i+1) * piecesize, size)
                sseqstr = seqstr[rbeg:rend]
                sid = "%s-%d-%d" % (id, bbeg+rbeg, bbeg+rend-1)
                rcd = SeqRecord(Seq(sseqstr), id = sid, description = '')
                rcds.append(rcd)
                #print "      " + sid
            SeqIO.write(rcds, fho, "fasta")
        else:
            SeqIO.write(seq, fho, "fasta")
    fhi.close()
    fho.close()

def merge(args):
    cfg = args.cfg
    for line in must_open(cfg):
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        (org, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            print os.access(fi, os.F_OK)
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
    fhc.close()
    return (orgs, fis)
def merge_seqs(fis, fids, fo):
    print "  merging input files to %s" % fo
    seqs = []
    for i in range(0,len(fids)):
        handle = 0
        if (fis[i].endswith(".gz")):
            handle = gzip.open(fis[i], "rb")
        else:
            handle = open(fis[i], "rU")
        seq_it = SeqIO.parse(handle, "fasta")
        handle.close

        seqs1 = [SeqRecord(rcd.seq, id = fids[i] + "|" + rcd.id,
            description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()


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
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    alns = AlignIO.parse(fhi, "fasta")
    AlignIO.write(alns, fho, "clustal")
    fhi.close()
    fho.close()

def rmdot(args):
    from string import maketrans
    tt = maketrans(".", "-")
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith('>'):
            fho.write(line)
        else:
            fho.write(line.translate(tt))
    fhi.close()
    fho.close()

def cleanid(args):
    fhi = open(args.fi, "r")
    fho = open(args.fo, "w")
    for line in fhi:
        if line.startswith(">"):
            fho.write(line.rstrip(":.\n")+"\n")
        else:
            fho.write(line)
    fhi.close()
    fho.close()

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

    sp3 = sp.add_parser("splitlong",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'break long fasta record into small pieces'
    )
    sp3.add_argument('fi', help = 'input file (fasta)')
    sp3.add_argument('--mode', type = int, default = 1, choices = [1, 2],
            help = 'split mode: 1 [100kb chunks], 2 [240 pieces]'
    )
    sp3.set_defaults(func = splitlong)

    sp3 = sp.add_parser("merge", help = 'merge multiple fasta files and update IDs')
    sp3.add_argument('cfg', help = 'config file (a text file with identifier followed by the absolute path of fasta in each line)'
    sp3.set_defaults(func = merge)
 
    sp9 = sp.add_parser("gaps",
            help = "report gap ('N's) locations in fasta sequences",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp9.add_argument('fi', help = 'input file (fasta)')
    sp9.add_argument('--gap', type = int, default = 10, help = 'min gap size')
    sp9.set_defaults(func = gaps)

    sp11 = sp.add_parser("fas2aln", help = 'convert fasta alignment file to clustal format')
    sp11.add_argument('fi', help = 'input alignment (.fas)')
    sp11.add_argument('fo', help = 'output alignment (.aln)')
    sp11.set_defaults(func = fas2aln)

    sp31 = sp.add_parser("rmdot", help = 'replace periods (.) in an alignment fasta by dashes (-)')
    sp31.add_argument('fi', help = 'input fasta file')
    sp31.add_argument('fo', help = 'output fasta file')
    sp31.set_defaults(func = rmdot)
    
    sp32 = sp.add_parser("cleanid", help = 'clean sequence IDs in a fasta file')
    sp32.add_argument('fi', help = 'input fasta file')
    sp32.add_argument('fo', help = 'output fasta file')
    sp32.set_defaults(func = cleanid)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Processing fastq files
"""

import os.path as op
import sys
import re
import logging
import json
import gzip

from itertools import islice

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from maize.formats.base import must_open, DictFile
from maize.utils.cbook import percentage
from maize.apps.base import sh, which, mkdir, need_update

FastqExt = ("fastq", "fq")
qual_offset = lambda x: 33 if x == "sanger" else 64
allowed_dialect_conversions = {
    ">=1.8": ("<1.8"),
    "sra": ("<1.8"),
}

class FastqLite (object):
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))

    def rc(self):
        self.seq = rc(self.seq)
        self.qual = self.qual[::-1]

class FastqRecord (object):
    def __init__(self, fh, offset=0, key=None):
        self.name = self.header = fh.readline()
        if not self.name:
            return
        self.name = self.name.split()[0]
        self.seq = fh.readline().rstrip()
        self.l3 = fh.readline().rstrip()
        self.qual = fh.readline().rstrip()
        if offset != 0:
            self.qual = "".join(chr(ord(x) + offset) for x in self.qual)
        self.length = len(self.seq)
        assert self.length == len(self.qual), \
                "length mismatch: seq(%s) and qual(%s)" % (self.seq, self.qual)
        if key:
            self.name = key(self.name)

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))

    def __len__(self):
        return self.length

    @property
    def quality(self):
        return [ord(x) for x in self.qual]

class FastqHeader(object):

    def __init__(self, row):
        header = row.strip().split(" ")
        self.readId, self.readLen, self.readNum = None, None, None
        self.multiplexId = 0
        self.paired = False
        if len(header) == 3 and "length" in header[2]:
            self.dialect = "sra"
            self.readId = header[0].lstrip('@')
            m = re.search("length\=(\d+)", header[2])
            if m:
                self.readLen = m.group(1)
            h = header[1].split(":")

            self.instrument = h[0]
            if len(h) == 7:
                self.runId, self.flowcellId = int(h[1]), h[2]
                self.laneNum, self.tileNum = int(h[3]), int(h[4])
                self.xPos, self.yPos = h[5], h[6]
            else:
                self.runId, self.flowcellId = None, None
                self.laneNum, self.tileNum = int(h[1]), int(h[2])
                self.xPos, self.yPos = h[3], h[4]
        else:
            h = header[0].split(":")
            self.instrument = h[0].lstrip('@')
            if len(header) == 2 and header[1].find(":"):
                self.dialect = ">=1.8"  # Illumina Casava 1.8+ format

                self.runId = int(h[1])
                self.flowcellId = h[2]
                self.laneNum = int(h[3])
                self.tileNum = int(h[4])
                self.xPos = int(h[5])
                self.yPos = h[6]
                if re.search("/", self.yPos):
                    self.paired = True
                    self.yPos, self.readNum = self.yPos.split("/")

                a = header[1].split(":")
                self.readNum = int(a[0])
                self.isFiltered = a[1]
                self.controlNum = int(a[2])
                self.barcode = a[3]
            else:
                self.dialect = "<1.8"   # Old Illumina Casava format (< 1.8)
                self.laneNum = int(h[1])
                self.tileNum = int(h[2])
                self.xPos = int(h[3])
                self.yPos = h[4]
                m = re.search(r"(\d+)(#\S+)\/(\d+)", self.yPos)
                if m:
                    self.paired = True
                    self.yPos, self.multiplexId, self.readNum = \
                            m.group(1), m.group(2), m.group(3)


    def __str__(self):
        if self.dialect == "sra":
            h0 = self.readId
            if self.readNum:
                h0 += "/{0}".format(self.readNum)

            h1elems = [self.instrument, self.laneNum, self.tileNum, \
                    self.xPos, self.yPos]
            if self.runId and self.flowcellId:
                h1elems[1:1] = [self.runId, self.flowcellId]
            h1 = ":".join(str(x) for x in h1elems)
            h2 = "length={0}".format(self.readLen)

            return "@{0} {1} {2}".format(h0, h1, h2)
        elif self.dialect == ">=1.8":
            yPos = "{0}/{1}".format(self.yPos, self.readNum) if self.paired \
                    else self.yPos

            h0 = ":".join(str(x) for x in (self.instrument, self.runId, \
                    self.flowcellId, self.laneNum, self.tileNum, \
                    self.xPos, yPos))
            h1 = ":".join(str(x) for x in (self.readNum, self.isFiltered, \
                    self.controlNum, self.barcode))

            return "@{0} {1}".format(h0, h1)
        else:
            yPos = "{0}#{1}/{2}".format(self.yPos, self.multiplexId, \
                    self.readNum) if self.paired else self.yPos
            h0 = ":".join(str(x) for x in (self.instrument, self.laneNum, \
                    self.tileNum, self.xPos, yPos))

            return "@{0}".format(h0)

    def format_header(self, dialect=None, tag=None):
        if dialect:
            if self.dialect == dialect:
                logging.error("Input and output dialect are the same")
            elif dialect not in allowed_dialect_conversions[self.dialect]:
                logging.error("Error: Cannot convert from `{0}` to `{1}` dialect".format(self.dialect, dialect))
                logging.error("Allowed conversions: {0}".format(json.dumps(allowed_dialect_conversions, indent=4)))
                sys.exit()
            else:
                self.dialect = dialect

        if tag:
            readNum = tag.split("/")[1]
            self.readNum = readNum
            self.paired = True

        return str(self)

def pairspf(pp, commonprefix=True):
    if commonprefix:
        pf = op.commonprefix(pp).rstrip("._-")
    else:
        pf = min(pp)
    pf = op.basename(pf)
    if not pf.strip():
        pf = op.basename(pp[0])
    return pf

def iter_fastq(filename, offset=0, key=None):
    if isinstance(filename, str):
        logging.debug("Read file `{0}`".format(filename))
        fh = must_open(filename)
    else:
        fh = filename

    while True:
        rec = FastqRecord(fh, offset=offset, key=key)
        if not rec.name:
            break
        yield rec
    yield None  # sentinel

def uniq(args):
    """
    %prog uniq fastqfile

    Retain only first instance of duplicate reads. Duplicate is defined as
    having the same read name.
    """
    p = OptionParser(uniq.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    fw = must_open(args.outfile, "w")
    nduplicates = nreads = 0
    seen = set()
    for rec in iter_fastq(fastqfile):
        nreads += 1
        if rec is None:
            break
        name = rec.name
        if name in seen:
            nduplicates += 1
            continue
        seen.add(name)
        print >> fw, rec
    logging.debug("Removed duplicate reads: {}".\
                  format(percentage(nduplicates, nreads)))

def suffix(args):
    """
    %prog suffix fastqfile CAG

    Filter reads based on suffix.
    """
    p = OptionParser(suffix.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastqfile, sf = args
    fw = must_open(args.outfile, "w")
    nreads = nselected = 0
    for rec in iter_fastq(fastqfile):
        nreads += 1
        if rec is None:
            break
        if rec.seq.endswith(sf):
            print >> fw, rec
            nselected += 1
    logging.debug("Selected reads with suffix {0}: {1}".\
                  format(sf, percentage(nselected, nreads)))

def calc_readlen(f, first):
    from maize.utils.cbook import SummaryStats

    L = []
    ai = iter_fastq(f)
    rec = ai.next()
    while rec:
        L.append(rec.length)
        if len(L) > first:
            break
        rec = ai.next()
    s = SummaryStats(L)

    return s

def is_fastq(f):
    fq = f.replace(".gz", "") if f.endswith(".gz") else f
    return fq.endswith((".fastq", ".fq"))

def readlen(args):
    """
    %prog readlen fastqfile

    Calculate read length, will only try the first N reads. Output min, max, and
    avg for each file.
    """
    p = OptionParser(readlen.__doc__)
    p.set_firstN()
    sp1.add_argument("--silent", default=False, action="store_true",
                 help="Do not print read length stats")
    sp1.add_argument("--nocheck", default=False, action="store_true",
                 help="Do not check file type suffix")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    f, = args
    if (not args.nocheck) and (not is_fastq(f)):
        logging.debug("File `{}` does not endswith .fastq or .fq".format(f))
        return 0

    s = calc_readlen(f, args.firstN)
    if not args.silent:
        print("\t".join(str(x) for x in (f, s.min, s.max, s.mean, s.median)))

    return int(s.max)

def fasta(args):
    """
    %prog fasta fastqfiles

    Convert fastq to fasta and qual file.
    """
    p = OptionParser(fasta.__doc__)
    sp1.add_argument("--seqtk", default=False, action="store_true",
                 help="Use seqtk to convert")
    p.set_outdir()
    p.set_outfile(outfile=None)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastqfiles = args
    outdir = args.outdir
    if outdir and outdir != ".":
        mkdir(outdir)

    fastqfile = fastqfiles[0]
    pf = op.basename(fastqfile)
    gzinput = pf.endswith(".gz")
    if gzinput:
        pf = pf.rsplit(".", 1)[0]

    pf, sf = pf.rsplit(".", 1)
    if sf not in ("fq", "fastq"):
        logging.debug("Assumed FASTA: suffix not `fq` or `fastq`")
        return fastqfile, None

    fastafile, qualfile = pf + ".fasta", pf + ".qual"
    outfile = args.outfile or fastafile
    outfile = op.join(outdir, outfile)
    if args.seqtk:
        if need_update(fastqfiles, outfile):
            for i, fastqfile in enumerate(fastqfiles):
                cmd = "seqtk seq -A {0} -L 30 -l 70".format(fastqfile)
                # First one creates file, following ones append to it
                sh(cmd, outfile=outfile, append=i)
        else:
            logging.debug("Outfile `{0}` already exists.".format(outfile))
        return outfile, None

    for fastqfile in fastqfiles:
        SeqIO.convert(fastqfile, "fastq", fastafile, "fasta")
        SeqIO.convert(fastqfile, "fastq", qualfile, "qual")

    return fastafile, qualfile

def first(args):
    """
    %prog first N fastqfile(s)

    Get first N reads from file.
    """
    from maize.apps.base import need_update

    p = OptionParser(first.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    N = int(args[0])
    nlines = N * 4
    fastqfiles = args[1:]
    fastqfile = fastqfiles[0]
    outfile = args.outfile
    if not need_update(fastqfiles, outfile):
        logging.debug("File `{0}` exists. Will not overwrite.".format(outfile))
        return

    gz = fastqfile.endswith(".gz")
    for fastqfile in fastqfiles:
        if gz:
            cmd = "zcat {0} | head -n {1}".format(fastqfile, nlines)
        else:
            cmd = "head -n {0} {1}".format(nlines, fastqfile)

        sh(cmd, outfile=args.outfile, append=True)

def FastqPairedIterator(read1, read2):
    if read1 == read2:
        p1fp = p2fp = must_open(read1)
    else:
        p1fp = must_open(read1)
        p2fp = must_open(read2)

    return p1fp, p2fp

def isHighQv(qs, qvchar, pct=90):
    cutoff = len(qs) * pct / 100
    highs = sum(1 for x in qs if x >= qvchar)
    return highs >= cutoff

def filter(args):
    """
    %prog filter paired.fastq

    Filter to get high qv reads. Use interleaved format (one file) or paired
    format (two files) to filter on paired reads.
    """
    p = OptionParser(filter.__doc__)
    sp1.add_argument("-q", dest="qv", default=20, type="int",
                 help="Minimum quality score to keep [default: %default]")
    sp1.add_argument("-p", dest="pct", default=95, type="int",
                 help="Minimum percent of bases that have [-q] quality "\
                 "[default: %default]")

    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    if len(args) == 1:
        r1 = r2 = args[0]
    else:
        r1, r2 = args

    qv = args.qv
    pct = args.pct

    offset = guessoffset([r1])
    qvchar = chr(offset + qv)
    logging.debug("Call base qv >= {0} as good.".format(qvchar))
    outfile = r1.rsplit(".", 1)[0] + ".q{0}.paired.fastq".format(qv)
    fw = open(outfile, "w")

    p1fp, p2fp = FastqPairedIterator(r1, r2)
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break

        b = list(islice(p2fp, 4))
        q1 = a[-1].rstrip()
        q2 = b[-1].rstrip()

        if isHighQv(q1, qvchar, pct=pct) and isHighQv(q2, qvchar, pct=pct):
            fw.writelines(a)
            fw.writelines(b)

def checkShuffleSizes(p1, p2, pairsfastq, extra=0):
    from maize.apps.base import getfilesize

    pairssize = getfilesize(pairsfastq)
    p1size = getfilesize(p1)
    p2size = getfilesize(p2)
    assert pairssize == p1size + p2size + extra, \
          "The sizes do not add up: {0} + {1} + {2} != {3}".\
          format(p1size, p2size, extra, pairssize)

def shuffle(args):
    """
    %prog shuffle p1.fastq p2.fastq

    Shuffle pairs into interleaved format.
    """
    p = OptionParser(shuffle.__doc__)
    p.set_tag()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    p1, p2 = args
    pairsfastq = pairspf((p1, p2)) + ".fastq"
    tag = args.tag

    p1fp = must_open(p1)
    p2fp = must_open(p2)
    pairsfw = must_open(pairsfastq, "w")
    nreads = 0
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break

        b = list(islice(p2fp, 4))
        if tag:
            name = a[0].rstrip()
            a[0] = name + "/1\n"
            b[0] = name + "/2\n"

        pairsfw.writelines(a)
        pairsfw.writelines(b)
        nreads += 2

    pairsfw.close()
    extra = nreads * 2 if tag else 0
    checkShuffleSizes(p1, p2, pairsfastq, extra=extra)

    logging.debug("File `{0}` verified after writing {1} reads.".\
                     format(pairsfastq, nreads))
    return pairsfastq

def split(args):
    """
    %prog split pairs.fastq

    Split shuffled pairs into `.1.fastq` and `.2.fastq`, using `sed`. Can work
    on gzipped file.

    <http://seqanswers.com/forums/showthread.php?t=13776>
    """
    from maize.apps.grid import Jobs

    p = OptionParser(split.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pairsfastq, = args
    gz = pairsfastq.endswith(".gz")
    pf = pairsfastq.replace(".gz", "").rsplit(".", 1)[0]
    p1 = pf + ".1.fastq"
    p2 = pf + ".2.fastq"

    cmd = "zcat" if gz else "cat"
    p1cmd = cmd + " {0} | sed -ne '1~8{{N;N;N;p}}'".format(pairsfastq)
    p2cmd = cmd + " {0} | sed -ne '5~8{{N;N;N;p}}'".format(pairsfastq)

    if gz:
        p1cmd += " | gzip"
        p2cmd += " | gzip"
        p1 += ".gz"
        p2 += ".gz"

    p1cmd += " > " + p1
    p2cmd += " > " + p2

    args = [(p1cmd, ), (p2cmd, )]
    m = Jobs(target=sh, args=args)
    m.run()

    checkShuffleSizes(p1, p2, pairsfastq)

def guessoffset(args):
    """
    %prog guessoffset fastqfile

    Guess the quality offset of the fastqfile, whether 33 or 64.
    See encoding schemes: <http://en.wikipedia.org/wiki/FASTQ_format>

      SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS...............................
      ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
      LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL...............................
      !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
      |                         |    |        |                              |
     33                        59   64       73                            104

     S - Sanger        Phred+33,  raw reads typically (0, 40)
     X - Solexa        Solexa+64, raw reads typically (-5, 40)
     I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
     J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     L - Illumina 1.8+ Phred+33,  raw reads typically (0, 40)
        with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
    """
    p = OptionParser(guessoffset.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    ai = iter_fastq(fastqfile)
    rec = ai.next()
    offset = 64
    while rec:
        quality = rec.quality
        lowcounts = len([x for x in quality if x < 59])
        highcounts = len([x for x in quality if x > 74])
        diff = highcounts - lowcounts
        if diff > 10:
            break
        elif diff < -10:
            offset = 33
            break
        rec = ai.next()

    if offset == 33:
        print >> sys.stderr, "Sanger encoding (offset=33)"
    elif offset == 64:
        print >> sys.stderr, "Illumina encoding (offset=64)"

    return offset

def format(args):
    """
    %prog format fastqfile

    Format FASTQ file. Currently provides option to convert FASTQ header from
    one dialect to another.
    """
    p = OptionParser(format.__doc__)

    sp1.add_argument("--convert", default=None, choices=[">=1.8", "<1.8", "sra"],
                help="Convert fastq header to a different format" +
                " [default: %default]")
    p.set_tag(specify_tag=True)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    ai = iter_fastq(fastqfile)
    rec = ai.next()
    dialect = None
    while rec:
        h = FastqHeader(rec.header)
        if not dialect:
            dialect = h.dialect
            logging.debug("Input fastq dialect: `{0}`".format(dialect))
            if args.convert:
                logging.debug("Output fastq dialect: `{0}`".format(args.convert))

        rec.name = h.format_header(dialect=args.convert, tag=args.tag)

        print(rec)
        rec = ai.next()

def some(args):
    """
    %prog some idsfile afastq [bfastq]

    Select a subset of the reads with ids present in the idsfile.
    `bfastq` is optional (only if reads are paired)
    """
    p = OptionParser(some.__doc__)
    opts, args = p.parse_args(args)

    if len(args) not in (2, 3):
        sys.exit(not p.print_help())

    idsfile, afastq, = args[:2]
    bfastq = args[2] if len(args) == 3 else None

    ids = DictFile(idsfile, valuepos=None)

    ai = iter_fastq(open(afastq))
    arec = ai.next()
    if bfastq:
        bi = iter_fastq(open(bfastq))
        brec = bi.next()

    while arec:
        if arec.name[1:] in ids:
            print(arec)
            if bfastq:
                print(brec)

        arec = ai.next()
        if bfastq:
            brec = bi.next()

def trim(args):
    """
    %prog trim fastqfile

    Wraps `fastx_trimmer` to trim from begin or end of reads.
    """
    p = OptionParser(trim.__doc__)
    sp1.add_argument("-f", dest="first", default=0, type="int",
            help="First base to keep. Default is 1.")
    sp1.add_argument("-l", dest="last", default=0, type="int",
            help="Last base to keep. Default is entire read.")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    obfastqfile = op.basename(fastqfile)
    fq = obfastqfile.rsplit(".", 1)[0] + ".ntrimmed.fastq"
    if fastqfile.endswith(".gz"):
        fq = obfastqfile.rsplit(".", 2)[0] + ".ntrimmed.fastq.gz"

    cmd = "fastx_trimmer -Q33 "
    if args.first:
        cmd += "-f {0.first} ".format(opts)
    if args.last:
        cmd += "-l {0.last} ".format(opts)

    sh(cmd, infile=fastqfile, outfile=fq)

def catread(args):
    """
    %prog catread fastqfile1 fastqfile2

    Concatenate paired end reads into one. Useful for example to do single-end
    mapping and perform filtering on the whole read pair level.
    """
    p = OptionParser(catread.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    r1, r2 = args
    p1fp, p2fp = FastqPairedIterator(r1, r2)
    outfile = pairspf((r1, r2)) + ".cat.fastq"
    fw = must_open(outfile, "w")
    while True:
        a = list(islice(p1fp, 4))
        if not a:
            break
        atitle, aseq, _, aqual = a
        btitle, bseq, _, bqual = list(islice(p2fp, 4))
        print >> fw, "\n".join((atitle.strip(), aseq.strip() + bseq.strip(), \
                                "+", aqual.strip() + bqual.strip()))

def splitread(args):
    """
    %prog splitread fastqfile

    Split fastqfile into two read fastqfiles, cut in the middle.
    """
    p = OptionParser(splitread.__doc__)
    sp1.add_argument("-n", dest="n", default=76, type="int",
            help="Split at N-th base position [default: %default]")
    sp1.add_argument("--rc", default=False, action="store_true",
            help="Reverse complement second read [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pairsfastq, = args

    base = op.basename(pairsfastq).split(".")[0]
    fq1 = base + ".1.fastq"
    fq2 = base + ".2.fastq"
    fw1 = must_open(fq1, "w")
    fw2 = must_open(fq2, "w")

    fp = must_open(pairsfastq)
    n = args.n
    minsize = n * 8 / 5

    for name, seq, qual in FastqGeneralIterator(fp):
        if len(seq) < minsize:
            logging.error("Skipping read {0}, length={1}".format(name, len(seq)))
            continue

        name = "@" + name
        rec1 = FastqLite(name, seq[:n], qual[:n])
        rec2 = FastqLite(name, seq[n:], qual[n:])
        if args.rc:
            rec2.rc()

        print >> fw1, rec1
        print >> fw2, rec2

    logging.debug("Reads split into `{0},{1}`".format(fq1, fq2))
    fw1.close()
    fw2.close()

def size(args):
    """
    %prog size fastqfile

    Find the total base pairs in a list of fastq files
    """
    p = OptionParser(size.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    total_size = total_numrecords = 0
    for f in args:
        cur_size = cur_numrecords = 0
        for rec in iter_fastq(f):
            if not rec:
                break
            cur_numrecords += 1
            cur_size += len(rec)

        print(" ".join(str(x) for x in (op.basename(f), cur_numrecords, cur_size)))
        total_numrecords += cur_numrecords
        total_size += cur_size

    if len(args) > 1:
        print(" ".join(str(x) for x in ("Total", total_numrecords, total_size)))

def convert(args):
    """
    %prog convert in.fastq

    illumina fastq quality encoding uses offset 64, and sanger uses 33. This
    script creates a new file with the correct encoding. Output gzipped file if
    input is also gzipped.
    """
    p = OptionParser(convert.__doc__)
    p.set_phred()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    infastq, = args
    phred = args.phred or str(guessoffset([infastq]))
    ophred = {"64": "33", "33": "64"}[phred]

    gz = infastq.endswith(".gz")
    outfastq = infastq.rsplit(".", 1)[0] if gz else infastq
    pf, sf = outfastq.rsplit(".", 1)
    outfastq = "{0}.q{1}.{2}".format(pf, ophred, sf)
    if gz:
        outfastq += ".gz"

    fin = "illumina" if phred == "64" else "sanger"
    fout = "sanger" if phred == "64" else "illumina"

    seqret = "seqret"
    if infastq.endswith(".gz"):
        cmd = "zcat {0} | ".format(infastq)
        cmd += seqret + " fastq-{0}::stdin fastq-{1}::stdout".\
                format(fin, fout)
    else:
        cmd = seqret + " fastq-{0}::{1} fastq-{2}::stdout".\
                format(fin, infastq, fout)

    sh(cmd, outfile=outfastq)

    return outfastq

def pairinplace(args):
    """
    %prog pairinplace bulk.fastq

    Pair up the records in bulk.fastq by comparing the names for adjancent
    records. If they match, print to bulk.pairs.fastq, else print to
    bulk.frags.fastq.
    """
    from maize.utils.iter import pairwise

    p = OptionParser(pairinplace.__doc__)
    p.set_rclip()
    p.set_tag()
    sp1.add_argument("--base",
                help="Base name for the output files [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    base = args.base or op.basename(fastqfile).split(".")[0]

    frags = base + ".frags.fastq"
    pairs = base + ".pairs.fastq"
    if fastqfile.endswith(".gz"):
        frags += ".gz"
        pairs += ".gz"

    fragsfw = must_open(frags, "w")
    pairsfw = must_open(pairs, "w")

    N = args.rclip
    tag = args.tag
    strip_name = (lambda x: x[:-N]) if N else None

    fh_iter = iter_fastq(fastqfile, key=strip_name)
    skipflag = False  # controls the iterator skip
    for a, b in pairwise(fh_iter):
        if b is None:  # hit the eof
            break

        if skipflag:
            skipflag = False
            continue

        if a.name == b.name:
            if tag:
                a.name += "/1"
                b.name += "/2"
            print >> pairsfw, a
            print >> pairsfw, b
            skipflag = True
        else:
            print >> fragsfw, a

    # don't forget the last one, when b is None
    if not skipflag:
        print >> fragsfw, a

    logging.debug("Reads paired into `%s` and `%s`" % (pairs, frags))
    return pairs

def fromsra(args):
    """
    %prog fromsra srafile

    Convert sra file to fastq using the sratoolkit `fastq-dump`
    """
    p = OptionParser(fromsra.__doc__)
    sp1.add_argument("--paired", default=False, action="store_true",
            help="Specify if library layout is paired-end " + \
                 "[default: %default]")
    sp1.add_argument("--compress", default=None, choices=["gzip", "bzip2"],
            help="Compress output fastq files [default: %default]")
    p.set_outdir()
    p.set_grid()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    srafile, = args
    paired = args.paired
    compress = args.compress
    outdir = args.outdir

    script_path = which("fastq-dump")
    if not script_path:
        logging.error("Cannot find `fastq-dump` in the PATH")
        sys.exit()

    cmd = [script_path]
    if compress:
        cmd.append("--{0}".format(compress))
    if paired:
        cmd.append("--split-files")
    if outdir:
        cmd.append("--outdir {0}".format(outdir))
    cmd.append(srafile)

    outcmd = " ".join(cmd)
    sh(outcmd, grid=args.grid)

class Error(Exception):
    pass

class Line(str):
    """A line of text with associated filename and line number."""
    def error(self, message):
        """Return an error relating to this line."""
        return Error("{0}({1}): {2}\n{3}"
                     .format(self.filename, self.lineno, message, self))

class Lines(object):
    """Lines(filename, iterator) wraps 'iterator' so that it yields Line
    objects, with line numbers starting from 1. 'filename' is used in
    error messages.

    """
    def __init__(self, filename, iterator):
        self.filename = filename
        self.lines = enumerate(iterator, start=1)

    def __iter__(self):
        return self

    def __next__(self):
        lineno, s = next(self.lines)
        s = s.decode('utf-8')
        line = Line(s)
        line.filename = self.filename
        line.lineno = lineno
        return line

    # For compatibility with Python 2.
    next = __next__

def read_fastq(filename, iterator):
    """Read FASTQ data from 'iterator' (which may be a file object or any
    other iterator that yields strings) and generate tuples (sequence
    name, sequence data, quality data). 'filename' is used in error
    messages.

    """
    # This implementation follows the FASTQ specification given here:
    # <http://nar.oxfordjournals.org/content/38/6/1767.full>
    import re
    at_seqname_re = re.compile(r'@(.+)$')
    sequence_re = re.compile(r'[!-*,-~]*$')
    plus_seqname_re = re.compile(r'\+(.*)$')
    quality_re = re.compile(r'[!-~]*$')

    lines = Lines(filename, iterator)
    for line in lines:
        # First line of block is @<seqname>.
        m = at_seqname_re.match(line)
        if not m:
            raise line.error("Expected @<seqname> but found:")
        seqname = m.group(1)
        try:
            # One or more lines of sequence data.
            sequence = []
            for line in lines:
                m = sequence_re.match(line)
                if not m:
                    break
                sequence.append(m.group(0))
            if not sequence:
                raise line.error("Expected <sequence> but found:")

            # The line following the sequence data consists of a plus
            # sign and an optional sequence name (if supplied, it must
            # match the sequence name from the start of the block).
            m = plus_seqname_re.match(line)
            if not m:
                raise line.error("Expected +[<seqname>] but found:")
            if m.group(1) not in ['', seqname]:
                raise line.error("Expected +{} but found:".format(seqname))

            # One or more lines of quality data, containing the same
            # number of characters as the sequence data.
            quality = []
            n = sum(map(len, sequence))
            while n > 0:
                line = next(lines)
                m = quality_re.match(line)
                if not m:
                    raise line.error("Expected <quality> but found:")
                n -= len(m.group(0))
                if n < 0:
                    raise line.error("<quality> is longer than <sequence>:")
                quality.append(m.group(0))

            yield seqname, ''.join(sequence), ''.join(quality)

        except StopIteration:
            raise line.error("End of input before sequence was complete:")

def breakread(args):
    fhi = must_open(args.fi)
    fo1 = "%s_1.fq.gz" % fo
    fo2 = "%s_2.fq.gz" % fo
    fho1 = gzip.open(fo1, "wb")
    fho2 = gzip.open(fo2, "wb")

    for (seqid, seq, qual) in read_fastq(args.fi, fhi):
        assert len(seq) == readlen * 2 and len(qual) == readlen * 2, \
                "%s: seq[%d] qual[%d] not %d" % \
                (seqid, len(seq), len(qual), readlen)
        eles = seqid.split(" ")
        if len(eles) > 2: seqid = " ".join(eles[0:2])
        seq1, seq2 = seq[0:readlen], seq[readlen:readlen*2]
        qual1, qual2 = qual[0:readlen], qual[readlen:readlen*2]
        fho1.write(("@%s\n%s\n+\n%s\n" % (seqid, seq1, qual1)).encode('utf8'))
        fho2.write(("@%s\n%s\n+\n%s\n" % (seqid, seq2, qual2)).encode('utf8'))

def UMIcount(args):
    """
    %prog UMIcount fastqfile

    Report number of occurances of each unique UMI
    """
    fhi = must_open(args.fi)
    if args.fi.endswith(".gz"):
        fhi = gzip.open(args.fi, "r")

    ud = dict()
    for (seqid, seq, qual) in read_fastq(args.fi, fhi):
        umi = seqid.split(" ")[1].split("+")[1]
        if umi in ud:
            ud[umi] += 1
        else:
            ud[umi] = 1

    fho = must_open(args.fo, 'w')
    for umi, cnt in ud.items():
        fho.write("%s\t%s\n" % (umi, cnt))

    logging.debug("{} UMIs detected".format(len(ud)))

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fastq utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("break", help = "break each fastq seq into two seqs of equal lengths")
    sp1.add_argument('fi', help = 'input file (*.fastq or *.fastq.gz)')
    sp1.add_argument('fo', help = 'output prefix (*_1.fq.gz and *_2.fq.gz)')
    sp1.add_argument('readlen', type = int, help = 'read length')
    sp1.set_defaults(func = breakread)

    sp1 = sp.add_parser('UMIcount', help='count occurance of UMIs in file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input file (*.fastq or *.fastq.gz)')
    sp1.add_argument('fo', help = 'output table of UMI occurances (*.tsv)')
    sp1.set_defaults(func = UMIcount)

    sp1 = sp.add_parser('size', help='total base pairs in the fastq files',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = size)

    sp1 = sp.add_parser('shuffle', help='shuffle paired reads into the same file interleaved',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = shuffle)

    sp1 = sp.add_parser('split', help='split paired reads into two files',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = split)

    sp1 = sp.add_parser('splitread', help='split appended reads (from JGI)',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = splitread)
    
    sp1 = sp.add_parser('catread', help='cat pairs together (reverse of splitread)',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = catread)
    
    sp1 = sp.add_parser('pairinplace', help='collect pairs by checking adjacent ids',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = pairinplace)
    
    sp1 = sp.add_parser('convert', help='convert between illumina and sanger offset',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = convert)
    
    sp1 = sp.add_parser('first', help='get first N reads from file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = first)
    
    sp1 = sp.add_parser('filter', help='filter to get high qv reads',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = filter)
    
    sp1 = sp.add_parser('suffix', help='filter reads based on suffix',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = suffix)
    
    sp1 = sp.add_parser('trim', help='trim reads using fastx_trimmer',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = trim)
    
    sp1 = sp.add_parser('some', help='select a subset of fastq reads',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = some)
    
    sp1 = sp.add_parser('guessoffset', help='guess the quality offset of the fastq records',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = guessoffset)
    
    sp1 = sp.add_parser('readlen', help='calculate read length',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = readlen)
    
    sp1 = sp.add_parser('format', help='format fastq file, convert header from casava 1.8+ to older format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = format)
    
    sp1 = sp.add_parser('fasta', help='convert fastq to fasta and qual file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = fasta)
    
    sp1 = sp.add_parser('fromsra', help='convert sra to fastq using `fastq-dump`',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = fromsra)

    sp1 = sp.add_parser('uniq', help='retain only first instance of duplicate (by name) reads',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = uniq)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

if __name__ == '__main__':
    main()


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

FastqExt = ("fastq", "fq")

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



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'fastq utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("break", help = "break each fastq seq into two seqs of equal lengths")
    sp1.add_argument('fi', help = 'input file (*fastq)')
    sp1.add_argument('fo', help = 'output prefix (*_1.fq.gz and *_2.fq.gz)')
    sp1.add_argument('readlen', type = int, help = 'read length')
    sp1.set_defaults(func = breakread)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


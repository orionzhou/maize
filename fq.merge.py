#!/usr/bin/env python
import os
import os.path as op
import sys
import gzip

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

def read_fasta(filename, iterator):
    import re
    seqname_re = re.compile(r'>(.+)$')
    sequence_re = re.compile(r'[!-=?-~]*$')

    lines = Lines(filename, iterator)
    seqname, sequence = None, []
    for line in lines:
        m = seqname_re.match(line)
        if m:
            if seqname:
                if not sequence:
                    raise line.error("empty seq for %s" % seqname)
                yield seqname, ''.join(sequence)
            seqname, sequence = m.group(1), []
        else:
            m = sequence_re.match(line)
            if not m:
                raise line.error("Expected <sequence> but found:")
            sequence.append(m.group(0))
    if seqname:
        yield seqname, ''.join(sequence)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
            'fi1', help = 'input fasta 1'
    )
    parser.add_argument(
            'fi2', help = 'input fasta 2'
    )
    parser.add_argument(
            'fo', help = 'output fasta file'
    )
    parser.add_argument(
            '--nosuffix', action = "store_true", 
            help = 'disable adding *.1 and *.2 sufix'
    )
    parser.add_argument(
            '--join', action = "store_true", 
            help = 'join seqs of two reads to make one long read'
    )
    args = parser.parse_args()
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



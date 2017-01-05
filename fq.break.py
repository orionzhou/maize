#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
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
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'break each fastq seq into two seqs of equal lengths'
    )
    parser.add_argument(
        'fi', help = 'input file (*.fastq)'
    )
    parser.add_argument(
        'fo', help = 'output prefix (*_1.fastq.gz and *_2.fastq.gz)'
    )
    parser.add_argument(
        'readlen', type = int, help = 'read length'
    )
    args = parser.parse_args()
    (fi, fo, readlen) = (args.fi, args.fo, args.readlen)
    assert op.isfile(fi), "cannot read %s" % fi
    
    (root, ext) = op.splitext(fi)
    if ext == '.gz':
        fhi = gzip.open(fi, "rb")
    else:
        fhi = open(fi, "rb")

    fo1 = "%s_1.fastq.gz" % fo
    fo2 = "%s_2.fastq.gz" % fo
    fho1 = gzip.open(fo1, "wb")
    fho2 = gzip.open(fo2, "wb")
    
    for (seqid, seq, qual) in read_fastq(fi, fhi):
        assert len(seq) == readlen * 2 and len(qual) == readlen * 2, \
                "%s: seq[%d] qual[%d] not %d" % \
                (seqid, len(seq), len(qual), readlen)
        eles = seqid.split(" ")
        if len(eles) > 2: seqid = " ".join(eles[0:2])
        seq1, seq2 = seq[0:readlen], seq[readlen:readlen*2]
        qual1, qual2 = qual[0:readlen], qual[readlen:readlen*2]
        fho1.write(("@%s\n%s\n+\n%s\n" % (seqid, seq1, qual1)).encode('utf8'))
        fho2.write(("@%s\n%s\n+\n%s\n" % (seqid, seq2, qual2)).encode('utf8'))



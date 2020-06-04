#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import math
import sys
import logging

#from Bio import SeqIO
from itertools import groupby, islice, cycle

from maize.apps.base import eprint, sh, debug, mkdir
debug()

FastaExt = ("fasta", "fas", "fa", "fna", "cds", "pep", "faa", "fsa", "seq", "nt", "aa")
FastqExt = ("fastq", "fq")

class BaseFile (object):

    def __init__(self, filename):

        self.filename = filename
        if filename:
            logging.debug("Load file `{0}`".format(filename))

class LineFile (BaseFile, list):
    """
    Generic file parser for line-based files
    """
    def __init__(self, filename, comment=None, load=False):

        super(LineFile, self).__init__(filename)

        if load:
            fp = must_open(filename)
            self.lines = [l.strip() for l in fp if l[0]!=comment]
            logging.debug("Load {0} lines from `{1}`.".\
                        format(len(self.lines), filename))

class DictFile (BaseFile, dict):
    """
    Generic file parser for multi-column files, keyed by a particular index.
    """
    def __init__(self, filename, keypos=0, valuepos=1, delimiter=None,
                       strict=True, keycast=None, cast=None):

        super(DictFile, self).__init__(filename)
        self.keypos = keypos

        fp = must_open(filename)
        ncols = max(keypos, valuepos) + 1
        thiscols = 0
        for lineno, row in enumerate(fp):
            row = row.rstrip()
            atoms = row.split(delimiter)
            thiscols = len(atoms)
            if thiscols < ncols:
                action = "Aborted" if strict else "Skipped"

                msg = "Must contain >= {0} columns.  {1}.\n".format(ncols, action)
                msg += "  --> Line {0}: {1}".format(lineno + 1, row)
                logging.error(msg)
                if strict:
                    sys.exit(1)
                else:
                    continue

            key = atoms[keypos]
            value = atoms[valuepos] if (valuepos is not None) else atoms
            if keycast:
                key = keycast(key)
            if cast:
                value = cast(value)
            self[key] = value

        assert thiscols, "File empty"
        self.ncols = thiscols
        logging.debug("Imported {0} records from `{1}`.".\
                    format(len(self), filename))

class SetFile (BaseFile, set):

    def __init__(self, filename, column=-1, delimiter=None):
        super(SetFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if not row.strip():
                continue
            keys = [x.strip() for x in row.split(delimiter)]
            if column >= 0:
                keys = [keys[column]]
            self.update(keys)

class FileShredder (object):
    """
    Same as rm -f *
    """
    def __init__(self, filelist, verbose=True):

        filelist = [x for x in filelist if x and op.exists(x)]
        cmd = "rm -rf {0}".format(" ".join(filelist))
        sh(cmd, log=verbose)

class FileMerger (object):
    """
    Same as cat * > filename
    """
    def __init__(self, filelist, outfile):

        self.filelist = filelist
        self.outfile = outfile
        self.ingz = filelist[0].endswith(".gz")
        self.outgz = outfile.endswith(".gz")

    def merge(self, checkexists=False):
        outfile = self.outfile
        if checkexists and not need_update(self.filelist, outfile):
            logging.debug("File `{0}` exists. Merge skipped.".format(outfile))
            return

        files = " ".join(self.filelist)
        ingz, outgz = self.ingz, self.outgz
        if ingz and outgz:  # can merge gz files directly
            cmd = "cat {0} > {1}".format(files, outfile)
            sh(cmd)
        else:
            cmd = "zcat" if self.ingz else "cat"
            cmd += " " + files
            sh(cmd, outfile=outfile)

        return outfile

class FileSplitter (object):

    def __init__(self, filename, outputdir=None, format="fasta", mode="cycle"):
        self.filename = filename
        self.outputdir = outputdir
        self.mode = mode

        format = format or self._guess_format(filename)
        logging.debug("format is %s" % format)

        if format in ("fasta", "fastq"):
            self.klass = "seqio"
        elif format == "clust":
            self.klass = "clust"
        else:
            self.klass = "txt"

        self.format = format
        mkdir(outputdir)

    def _open(self, filename):

        if self.klass == "seqio":
            handle = SeqIO.parse(open(filename), self.format)
        elif self.klass == "clust":
            from jcvi.apps.uclust import ClustFile
            handle = iter(ClustFile(filename))
        else:
            handle = open(filename)
        return handle

    @property
    def num_records(self):
        handle = self._open(self.filename)
        return sum(1 for x in handle)

    def _guess_format(self, filename):
        root, ext = op.splitext(filename)
        ext = ext.strip(".")

        if ext in FastaExt:
            format = "fasta"
        elif ext in FastqExt:
            format = "fastq"
        else:
            format = "txt"
        return format

    def _batch_iterator(self, N=1):
        """Returns N lists of records.

        This can be used on any iterator, for example to batch up
        SeqRecord objects from Bio.SeqIO.parse(...), or to batch
        Alignment objects from Bio.AlignIO.parse(...), or simply
        lines from a file handle.

        This is a generator function, and it returns lists of the
        entries from the supplied iterator.  Each list will have
        batch_size entries, although the final list may be shorter.
        """
        batch_size = math.ceil(self.num_records / float(N))
        handle = self._open(self.filename)
        while True:
            batch = list(islice(handle, batch_size))
            if not batch:
                break
            yield batch

    @classmethod
    def get_names(cls, filename, N):
        root, ext = op.splitext(op.basename(filename))

        names = []
        pad0 = len(str(int(N - 1)))
        for i in xrange(N):
            name = "{0}_{1:0{2}d}{3}".format(root, i, pad0, ext)
            names.append(name)

        return names

    def write(self, fw, batch):
        if self.klass == "seqio":
            SeqIO.write(batch, fw, self.format)
        elif self.klass == "clust":
            for b in batch:
                print >> fw, b
        else:
            for line in batch:
                fw.write(line)
        return len(batch)

    def split(self, N, force=False):
        """
        There are two modes of splitting the records
        - batch: splitting is sequentially to records/N chunks
        - cycle: placing each record in the splitted files and cycles

        use `cycle` if the len of the record is not evenly distributed
        """
        mode = self.mode
        assert mode in ("batch", "cycle", "optimal")
        logging.debug("set split mode=%s" % mode)

        self.names = self.__class__.get_names(self.filename, N)
        if self.outputdir:
            self.names = [op.join(self.outputdir, x) for x in self.names]

        if not need_update(self.filename, self.names) and not force:
            logging.error("file %s already existed, skip file splitting" % \
                    self.names[0])
            return

        filehandles = [open(x, "w") for x in self.names]

        if mode == "batch":
            for batch, fw in zip(self._batch_iterator(N), filehandles):
                count = self.write(fw, batch)
                logging.debug("write %d records to %s" % (count, fw.name))

        elif mode == "cycle":
            handle = self._open(self.filename)
            for record, fw in izip(handle, cycle(filehandles)):
                count = self.write(fw, [record])

        elif mode == "optimal":
            """
            This mode is based on Longest Processing Time (LPT) algorithm:

            A simple, often-used algorithm is the LPT algorithm (Longest
            Processing Time) which sorts the jobs by its processing time and
            then assigns them to the machine with the earliest end time so far.
            This algorithm achieves an upper bound of 4/3 - 1/(3m) OPT.

            Citation: <http://en.wikipedia.org/wiki/Multiprocessor_scheduling>
            """
            endtime = [0] * N
            handle = self._open(self.filename)
            for record in handle:
                mt, mi = min((x, i) for (i, x) in enumerate(endtime))
                fw = filehandles[mi]
                count = self.write(fw, [record])
                endtime[mi] += len(record)

        for fw in filehandles:
            fw.close()

def longest_unique_prefix(query, targets, remove_self=True):
    """
    Find the longest unique prefix for filename, when compared against a list of
    filenames. Useful to simplify file names in a pool of files. See usage in
    formats.fasta.pool().
    """
    query = op.basename(query)
    targets = [op.basename(x) for x in targets]
    prefix_lengths = [len(op.commonprefix([query, name])) for name in targets]
    if remove_self and len(query) in prefix_lengths:
        prefix_lengths.remove(len(query))
    longest_length = max(prefix_lengths)
    return query[:longest_length + 1]

def check_exists(filename, oappend=False):
    """
    Avoid overwriting some files accidentally.
    """
    if op.exists(filename):
        if oappend:
            return oappend
        logging.error("`{0}` found, overwrite (Y/N)?".format(filename))
        overwrite = (raw_input() == 'Y')
    else:
        overwrite = True

    return overwrite

def timestamp():
    from datetime import datetime as dt
    return "{0}{1:02d}{2:02d}".format(dt.now().year, dt.now().month, dt.now().day)

def must_open(filename, mode="r", checkexists=False, skipcheck=False, \
            oappend=False):
    """
    Accepts filename and returns filehandle.

    Checks on multiple files, stdin/stdout/stderr, .gz or .bz2 file.
    """
    if isinstance(filename, list):
        assert "r" in mode

        if filename[0].endswith((".gz", ".bz2")):
            filename = " ".join(filename)  # allow opening multiple gz/bz2 files
        else:
            import fileinput
            return fileinput.input(filename)

    if filename.startswith("s3://"):
        from jcvi.utils.aws import pull_from_s3
        filename = pull_from_s3(filename)

    if filename in ("-", "stdin"):
        assert "r" in mode
        fp = sys.stdin

    elif filename == "stdout":
        assert "w" in mode
        fp = sys.stdout

    elif filename == "stderr":
        assert "w" in mode
        fp = sys.stderr

    elif filename == "tmp" and mode == "w":
        from tempfile import NamedTemporaryFile
        fp = NamedTemporaryFile(delete=False)

    elif filename.endswith(".gz"):
        if "r" in mode:
            cmd = "gunzip -c %s" % filename
            from subprocess import Popen, PIPE
            fp = Popen(cmd, bufsize = 1, stdout = PIPE, shell = True, universal_newlines = True).stdout
        elif "w" in mode:
            import gzip
            fp = gzip.open(filename, mode)

    elif filename.endswith(".bz2"):
        if "r" in mode:
            cmd = "bzcat -c %s" % filename
            from subprocess import Popen, PIPE
            fp = Popen(cmd, bufsize = 1, stdout = PIPE, shell = True, universal_newlines = True).stdout
        elif "w" in mode:
            import bz2
            fp = bz2.open(filename, mode)

    else:
        if checkexists:
            assert mode == "w"
            overwrite = (not op.exists(filename)) if skipcheck \
                        else check_exists(filename, oappend)
            if overwrite:
                if oappend:
                    fp = open(filename, "a")
                else:
                    fp = open(filename, "w")
            else:
                logging.debug("File `{0}` already exists. Skipped."\
                        .format(filename))
                return None
        else:
            fp = open(filename, mode)

    return fp

bash_shebang = "#!/bin/bash"
python_shebang = """#!/usr/bin/env python
# -*- coding: UTF-8 -*-"""

def write_file(filename, contents, meta=None, skipcheck=False, append=False, tee=False):
    if not meta:
        suffix = filename.rsplit(".", 1)[-1]
        if suffix == "sh":
            meta = "run script"
        elif suffix == "py":
            meta = "python script"
        else:
            meta = "file"

    meta_choices = ("file", "run script", "python script")
    assert meta in meta_choices, "meta must be one of {0}".\
                    format("|".join(meta_choices))

    contents = contents.strip()
    shebang = "\n"
    if "script" in meta:
        if not append:
            if meta == "run script":
                shebang = bash_shebang
            elif meta == "python script":
                shebang = python_shebang
        contents = "\n\n".join((shebang, contents))

    fw = must_open(filename, "w", checkexists=True, skipcheck=skipcheck, oappend=append)
    if fw:
        print >> fw, contents
        fw.close()
    if tee:
        print >> sys.stderr, contents

    fileop = "appended" if append else "written"
    message = "{0} {1} to `{2}`.".format(meta, fileop, filename)
    logging.debug(message.capitalize())
    if meta == "run script" and not append:
        sh("chmod u+x {0}".format(filename))

def read_until(handle, start):
    # read each line until a certain start, then puts the start tag back
    while 1:
        pos = handle.tell()
        line = handle.readline()
        if not line:
            break
        if line.startswith(start):
            handle.seek(pos)
            return

def read_block(handle, signal):
    """
    Useful for reading block-like file formats, for example FASTA or OBO file,
    such file usually startswith some signal, and in-between the signals are a
    record
    """
    signal_len = len(signal)
    it = (x[1] for x in groupby(handle,
        key=lambda row: row.strip()[:signal_len] == signal))
    found_signal = False
    for header in it:
        header = list(header)
        for h in header[:-1]:
            h = h.strip()
            if h[:signal_len] != signal:
                continue
            yield h, []  # Header only, no contents
        header = header[-1].strip()
        if header[:signal_len] != signal:
            continue
        found_signal = True
        seq = list(s.strip() for s in next(it))
        yield header, seq

    if not found_signal:
        handle.seek(0)
        seq = list(s.strip() for s in handle)
        yield None, seq

def is_number(s, cast=float):
    """
    Check if a string is a number. Use cast=int to check if s is an integer.
    """
    try:
        cast(s) # for int, long and float
    except ValueError:
        return False

    return True

def get_number(s, cast=int):
    """
    Try to get a number out of a string, and cast it.
    """
    import string
    d = "".join(x for x in str(s) if x in string.digits)
    return cast(d)

def flexible_cast(s):
    if is_number(s, cast=int):
        return int(s)
    elif is_number(s, cast=float):
        return float(s)
    return s

def ndigit(num):
    if num < 1:
        eprint("no digits: %g" % num)
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

def flatten(args):
    """
    %prog flatten filename > ids

    Convert a list of IDs (say, multiple IDs per line) and move them into one
    per line.

    For example, convert this, to this:
    A,B,C                    | A
    1                        | B
    a,4                      | C
                             | 1
                             | a
                             | 4

    If multi-column file with multiple elements per column, zip then flatten like so:
    A,B,C    2,10,gg         | A,2
    1,3      4               | B,10
                             | C,gg
                             | 1,4
                             | 3,na
    """
    from itertools import izip_longest

    tabfile, sep, zipsep = args.file, args.sep, args.zipsep

    fp = must_open(tabfile)
    for row in fp:
        if zipsep:
            row = row.rstrip()
            atoms = row.split(sep)
            frows = []
            for atom in atoms:
                frows.append(atom.split(zipsep))
            print("\n".join([zipsep.join(x) for x in list(izip_longest(*frows, fillvalue="na"))]))
        else:
            print(row.strip().replace(opts.sep, "\n"))

def split(args):
    """
    %prog split file outdir N

    Split file into N records. This allows splitting FASTA/FASTQ/TXT file
    properly at boundary of records. Split is useful for parallelization
    on input chunks.

    Option --mode is useful on how to break into chunks.
    1. chunk - chunk records sequentially, 1-100 in file 1, 101-200 in file 2, etc.
    2. cycle - chunk records in Round Robin fashion
    3. optimal - try to make split file of roughly similar sizes, using LPT
    algorithm. This is the default.
    """
    filename, outdir, N = args.file, args.outdir, args.N
    fs = FileSplitter(filename, outputdir=outdir,
                      format=args.format, mode=args.mode)

    if args.all:
        logging.debug("option -all override N")
        N = fs.num_records
    else:
        N = min(fs.num_records, int(N))
        assert N > 0, "N must be > 0"

    logging.debug("split file into %d chunks" % N)
    fs.split(N)

    return fs

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog = 'python -m maize.formats.base',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("split", 
            help = "split large file into N chunks",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp1.add_argument('file', help = 'input file')
    sp1.add_argument('outdir', help = 'output directory')
    sp1.add_argument('N', type = int, help = 'number of pieces')
    sp1.add_argument('--all', action = 'store_true', help = 'split all records')
    sp1.add_argument('--mode', choices = ['batch','cycle','optimal'], help = 'mode when splitting records')
    sp1.add_argument('--format', choices = ["fasta", "fastq", "txt", "clust"], help = 'input file format')
    sp1.set_defaults(func = split)

    sp2 = sp.add_parser("flatten", 
            help = "convert a list of IDs into one per line",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp2.add_argument('file', help = 'input file')
    sp2.add_argument('--sep', default = ',', help = 'input file separator')
    sp2.add_argument('--zipflatten', default = None, dest = 'zipsep',
            help = 'specify if columns should be zipped before flattening - if yes, specify output separator')

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

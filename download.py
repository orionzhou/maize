#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import wget

if sys.version_info[0] < 3:
    int_types = (int, long)
    urlopen = urllib.urlopen
else:
    int_types = (int,)
    basestring = str
    from urllib.request import urlopen

def download_file(src_ftp, dst_file, prt=sys.stdout, loading_bar=True):
    """Download specified file if necessary."""
    if os.path.isfile(dst_file):
        return
    do_gunzip = src_ftp[-3:] == '.gz' and dst_file[-3:] != '.gz'
    dst_wget = "{DST}.gz".format(DST=dst_file) if do_gunzip else dst_file
    # Write to stderr, not stdout so this message will be seen when running nosetests
    wget_msg = "wget.download({SRC} out={DST})\n".format(SRC=src_ftp, DST=dst_wget)
    sys.stderr.write("  {WGET}".format(WGET=wget_msg))
    if loading_bar:
        loading_bar = wget.bar_adaptive
    try:
        wget.download(src_ftp, out=dst_wget, bar=loading_bar)
        if do_gunzip:
            if prt is not None:
                prt.write("  gunzip {FILE}\n".format(FILE=dst_wget))
            gzip_open_to(dst_wget, dst_file)
    except IOError as errmsg:
        import traceback
        traceback.print_exc()
        sys.stderr.write("**FATAL cmd: {WGET}".format(WGET=wget_msg))
        sys.stderr.write("**FATAL msg: {ERR}".format(ERR=str(errmsg)))
        sys.exit(1)

def gzip_open_to(fin_gz, fout):
    """Unzip a file.gz file."""
    with gzip.open(fin_gz, 'rb') as zstrm:
        with  open(fout, 'wb') as ostrm:
            ostrm.write(zstrm.read())
    assert os.path.isfile(fout), "COULD NOT GUNZIP({G}) TO FILE({F})".format(G=fin_gz, F=fout)
    os.remove(fin_gz)

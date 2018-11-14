#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Manipulate PDF files, using PyPDF2 library.
"""

import os
import sys
import logging
import traceback

from PyPDF2 import PdfFileMerger, parse_filename_page_ranges
from PyPDF2.pagerange import PAGE_RANGE_HELP
from maize.formats.base import must_open
from maize.utils.natsort import natsorted

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = ''
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('cat', help='concatenate pages from pdf files into a single pdf file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('i', help = '')
    sp1.set_defaults(func = catx)
    
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

def cat(args):
    """
    %prog cat *.pdf -o output.pdf

    Concatenate pages from pdf files into a single pdf file.

    Page ranges refer to the previously-named file.
    A file not followed by a page range means all the pages of the file.

    PAGE RANGES are like Python slices.
            {page_range_help}
    EXAMPLES
        pdfcat -o output.pdf head.pdf content.pdf :6 7: tail.pdf -1
            Concatenate all of head.pdf, all but page seven of content.pdf,
            and the last page of tail.pdf, producing output.pdf.

        pdfcat chapter*.pdf >book.pdf
            You can specify the output file by redirection.

        pdfcat chapter?.pdf chapter10.pdf >book.pdf
            In case you don't want chapter 10 before chapter 2.
    """
    p = OptionParser(cat.__doc__.format(page_range_help=PAGE_RANGE_HELP))
    sp1.add_argument("--nosort", default=False, action="store_true",
                 help="Do not sort file names")
    sp1.add_argument("--cleanup", default=False, action="store_true",
                 help="Remove individual pdfs after merging")
    p.set_outfile()
    p.set_verbose(help="Show page ranges as they are being read")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    outfile = args.outfile
    if outfile in args:
        args.remove(outfile)

    if not args.nosort:
        args = natsorted(args)

    filename_page_ranges = parse_filename_page_ranges(args)
    verbose = args.verbose
    fw = must_open(outfile, "wb")

    merger = PdfFileMerger()
    in_fs = {}
    try:
        for (filename, page_range) in filename_page_ranges:
            if verbose:
                print >> sys.stderr, filename, page_range
            if filename not in in_fs:
                in_fs[filename] = open(filename, "rb")
            merger.append(in_fs[filename], pages=page_range)
    except:
        print >> sys.stderr, traceback.format_exc()
        print >> sys.stderr, "Error while reading " + filename
        sys.exit(1)
    merger.write(fw)
    fw.close()

    if args.cleanup:
        logging.debug("Cleaning up {} files".format(len(args)))
        for arg in args:
            os.remove(arg)

if __name__ == '__main__':
    main()

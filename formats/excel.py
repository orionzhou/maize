#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read and write EXCEL file.

Library dependency: pandas
"""

import os.path as op
import sys
import logging
import pandas as pd

from maize.apps.base import sh, mkdir

class ColorMatcher(object):

    def __init__(self):
        self.reset()

    def reset(self):
        self.unused_colors = set(self.xlwt_colors)
        # Never use black.
        self.unused_colors.discard((0, 0, 0))

    # Culled from a table at http://www.mvps.org/dmcritchie/excel/colors.htm
    xlwt_colors=[
        (0,0,0), (255,255,255), (255,0,0), (0,255,0), (0,0,255), (255,255,0),
        (255,0,255), (0,255,255), (0,0,0), (255,255,255), (255,0,0), (0,255,0),
        (0,0,255), (255,255,0), (255,0,255), (0,255,255), (128,0,0), (0,128,0),
        (0,0,128), (128,128,0), (128,0,128), (0,128,128), (192,192,192),
        (128,128,128), (153,153,255), (153,51,102), (255,255,204),
        (204,255,255), (102,0,102), (255,128,128), (0,102,204), (204,204,255),
        (0,0,128), (255,0,255), (255,255,0), (0,255,255), (128,0,128),
        (128,0,0), (0,128,128), (0,0,255), (0,204,255), (204,255,255),
        (204,255,204), (255,255,153), (153,204,255), (255,153,204),
        (204,153,255), (255,204,153), (51,102,255), (51,204,204), (153,204,0),
        (255,204,0), (255,153,0), (255,102,0), (102,102,153), (150,150,150),
        (0,51,102), (51,153,102), (0,51,0), (51,51,0), (153,51,0), (153,51,102),
        (51,51,153), (51,51,51)
    ]

    @staticmethod
    def color_distance(rgb1, rgb2):
        # Adapted from Colour metric by Thiadmer Riemersma,
        # http://www.compuphase.com/cmetric.htm
        rmean = (rgb1[0] + rgb2[0]) / 2
        r = rgb1[0] - rgb2[0]
        g = rgb1[1] - rgb2[1]
        b = rgb1[2] - rgb2[2]
        return (((512 + rmean) * r * r) / 256) + 4 * g * g \
            + (((767 - rmean) * b * b) / 256)

    def match_color_index(self, color):
        """Takes an "R,G,B" string or wx.Color and returns a matching xlwt
        color.
        """
        from maize.utils.webcolors import color_diff
        if isinstance(color, int):
            return color
        if color:
            if isinstance(color, basestring):
                rgb = map(int, color.split(','))
            else:
                rgb = color.Get()
            logging.disable(logging.DEBUG)
            distances = [color_diff(rgb, x) for x in self.xlwt_colors]
            logging.disable(logging.NOTSET)
            result = distances.index(min(distances))
            self.unused_colors.discard(self.xlwt_colors[result])
            return result

    def get_unused_color(self):
        """Returns an xlwt color index that has not been previously returned by
        this instance.  Attempts to maximize the distance between the color and
        all previously used colors.
        """
        if not self.unused_colors:
            # If we somehow run out of colors, reset the color matcher.
            self.reset()
        used_colors = [c for c in self.xlwt_colors if c not in self.unused_colors]
        result_color = max(self.unused_colors,
                           key=lambda c: min(self.color_distance(c, c2)
                                             for c2 in used_colors))
        result_index = self.xlwt_colors.index(result_color)
        self.unused_colors.discard(result_color)
        return result_index

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read and write excel files'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser('csv', help='Convert EXCEL to csv file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('excel', help = 'input excel file')
    sp1.add_argument('--sheet_name', default='Sheet1', help='worksheet name')
    sp1.add_argument('--sep', default=',', help='separator')
    sp1.set_defaults(func = csv)

    sp1 = sp.add_parser('tsv', help='Convert EXCEL to tsv file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('excel', help = 'input excel file')
    sp1.add_argument('--sheet_name', default='Sheet1', help='worksheet name')
    sp1.add_argument('--sep', default='\t', help='separator')
    sp1.set_defaults(func = csv)

    sp1 = sp.add_parser('tsvs', help='Convert all worksheets in EXCEL to tsv files',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('excel', help = 'input excel file')
    sp1.add_argument('--outdir', default='sheets', help='output directory')
    sp1.add_argument('--sep', default='\t', help='separator')
    sp1.set_defaults(func = tsvs)

    sp1 = sp.add_parser('fromcsv', help='Convert csv file to EXCEL',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('csv', help = 'input csv file')
    sp1.add_argument("--noheader", action="store_true",
                 help="Do not treat the first row as header")
    sp1.add_argument("--rgb", default=-1, type=int,
                 help="Show RGB color box")
    sp1.set_defaults(func = fromcsv)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

def fromcsv(args):
    """
    %prog fromcsv csvfile

    Convert csv file to EXCEL.
    """
    from csv import reader
    from xlwt import Workbook, easyxf
    from maize.formats.base import flexible_cast
    
    csvfile = args.csv
    header = not args.noheader
    rgb = args.rgb
    excelfile = csvfile.rsplit(".", 1)[0] + ".xls"

    data = []
    for row in reader(open(csvfile), delimiter=args.sep):
        data.append(row)

    w = Workbook()
    s = w.add_sheet(op.basename(csvfile))

    header_style = easyxf('font: bold on')
    if header:
        s.panes_frozen = True
        s.horz_split_pos = 1

    cm = ColorMatcher()
    for i, row in enumerate(data):
        for j, cell in enumerate(row):
            cell = flexible_cast(cell)
            if header and i == 0:
                s.write(i, j, cell, header_style)
            else:
                if j == rgb:
                    cix = cm.match_color_index(cell)
                    color_style = easyxf('font: color_index {0}'.format(cix))
                    s.write(i, j, cell, color_style)
                else:
                    s.write(i, j, cell)

    w.save(excelfile)
    logging.debug("File written to `{0}`.".format(excelfile))
    return excelfile

def csv(args):
    """
    %prog csv excelfile

    Convert EXCEL to csv file.
    """
    import pandas as pd

    excelfile = args.excel
    sep = args.sep
    sheet_name = args.sheet_name
    suf = '.tsv' if sep == '\t' else '.csv'
    csvfile = excelfile.rsplit(".", 1)[0] + suf
    df = pd.read_excel(excelfile, sheet_name=sheet_name, header=0)
    df.to_csv(csvfile, sep=sep, header=True, index=False)

def tsvs(args):
    """
    %prog tsvs excelfile

    Convert all worksheets in EXCEL to tsv files.
    """
    excelfile = args.excel
    odir = args.outdir
    sep = args.sep

    xl = pd.ExcelFile(excelfile)
    sheets = xl.sheet_names
    print("will convert %d sheets under %s" % (len(sheets), odir))
    mkdir(odir)

    suf = '.tsv' if sep == '\t' else '.csv'
    for sheet in sheets:
        fo = "%s/%s%s" % (odir, sheet, suf)
        print("    writing %s" % fo)
        df = pd.read_excel(excelfile, sheet_name=sheet, header=0)
        df.to_csv(fo, sep=sep, header=True, index=False)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import math
import sys
import logging

from jcvi.apps.base import sh, debug, mkdir

def ndigit(num):
    if num < 1:
        print("no digits: %g" % num)
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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog = 'python -m maize.formats.base',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

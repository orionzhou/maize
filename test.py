#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path
import argparse

def main(**kwargs):
    for key, value in kwargs.iteritems():
        print key, value
    cmd = '{0} {1}'.format(kwargs['program'], ' '.join(kwargs['infiles']))
    r = envoy.run(cmd)
    print r.std_out
    print r.std_err

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='JJJ')
    parser.add_argument('--version', action="version", version="%(prog)s 1.0")
    parser.add_argument('program', type=str, help='progname')
    parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, help='input')
    parser.add_argument('output', type=argparse.FileType('w'), default=sys.stdout, help='output')
    args = parser.parse_args()
    main(**vars(args))
 

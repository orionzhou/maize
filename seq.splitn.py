#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse
import numpy as np
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def digitof_number(num):
    if num < 1:
        print "no digits: %g", num
        sys.exit(1)
    digit = 0 
    while num >= 1:
        num /= 10.0
        digit += 1
    return digit
def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'run pyfasta to partition a set of fasta records evenly'
    )
    parser.add_argument(
        'fi', help = 'input file (fasta)'
    )
    parser.add_argument(
        'dirw', help = 'output directory'
    )
    parser.add_argument(
        '--pieces', dest = "n", default = 24, type = int, help = '# pieces to partition [24]'
    )
    nproc = int(os.environ['nproc'])
    parser.add_argument(
        '--cpu', dest = "cpu", default = nproc, type = int, help = 'number processors to use (default: all/%d)' % nproc
    )
    args = parser.parse_args()
    (fi, dirw) = (args.fi, args.dirw)
    (ncpu, n) = (args.cpu, args.n)
    (fi, dirw) = (op.realpath(fi), op.realpath(dirw))
    if not op.exists(dirw): 
        os.makedirs(dirw)
    else:
        os.system("rm -rf %s/*" % dirw)
    
    cdir = os.path.dirname(os.path.realpath(__file__))
    cwd = os.getcwd()
    os.chdir(dirw)

    os.system("ln -sf %s part.fas" % fi)
    cmd = "pyfasta split -n %d part.fas" % n
    os.system(cmd)
    os.system("rm part.fas part.fas.*")

    digit = digitof_number(n)
    sizes = []
    for i in range(0,n):
        fmt = "part.%%0%dd.fas" % digit
        fp = fmt % i
        sizes.append(os.stat(fp).st_size)
    sizes.sort()
    print "size range: %s - %s" % (sizeof_fmt(sizes[0]), sizeof_fmt(sizes[n-1]))


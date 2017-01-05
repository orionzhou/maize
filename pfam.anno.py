#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import argparse
from colorama import init, Fore, Back, Style

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Run Hmmscan with Pfam-A.hmm to annotate a set of genes'
    )
    parser.add_argument(
        'fi', help = 'Gene file to be annotated (*.Gtb)'
    )
    parser.add_argument(
        'fs', help = 'Protein sequences of input genes'
    )
    parser.add_argument(
        'fo', help = 'Output file (*.Gtb)'
    )
    parser.add_argument(
        '--fr', dest='fr', help='RepeatMasker output file (*.tbl)'
    )
    parser.add_argument(
        '--nproc', dest='nproc', type=int, default=1, help='Number of processors to run hmmscan (default: 1)'
    )
    parser.add_argument(
        '--pfam', dest='pfam', default='/home/springer/zhoux379/data/db/pfam_30/Pfam-A.hmm', help='HMM file (default: /home/springer/zhoux379/data/db/pfam_30/Pfam-A.hmm)'
    )
    args = parser.parse_args()
    (fi, fs, fo, fr, nproc, pfam) = (args.fi, args.fs, args.fo, args.fr, args.nproc, args.pfam)
    
    cmds = []
    (pre, ext) = op.splitext(fo)
    cmds.append("cd %s" % os.getcwd())
    cmds.append("hmmscan --cpu %d -o %s.1.txt %s %s" % (nproc, pre, pfam, fs))
    cmds.append("hmmc2htb.pl -i %s.1.txt -o %s.2.htb -m %s -s %s" % (pre, pre, pfam, fs))
    cmds.append("htb.qtile.pl -i %s.2.htb -o %s.3.htb" % (pre, pre))
    cmds.append("htb.filter.pl -i %s.3.htb -l 10 -e 0.01 -o %s.4.htb" % (pre, pre))
    cmds.append("cut -f2-4,6,7-9,11-13 %s.4.htb > %s.5.tsv" % (pre, pre))
    cmds.append("gtb.addpfam.pl -i %s -p %s.5.tsv -o %s.6.gtb" % (fi, pre, pre))
    cmds.append("gtb2bed.s.pl -i %s.6.gtb -o %s.7.bed" % (pre, pre))
    cmds.append("intersectBed -wao -a %s.7.bed -b %s > %s.8.bed" % (pre, fr, pre))
    cmds.append("gtb.addrm.pl -i %s.6.gtb -b %s.8.bed -o %s" % (pre, pre, fo))
    init()
    print(Fore.RED)
    print("Please copy the following lines into your job script:")
    print(Fore.GREEN + "\n".join(cmds))
    print(Style.RESET_ALL)




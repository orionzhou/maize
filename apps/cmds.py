#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from subprocess import Popen, PIPE, run
from jcvi.apps.base import sh, mkdir
from filecmp import cmp
from shutil import copy

def update_conda(args):
    envs1 = '''base snk work
        blast hisat2 bismark
        alfred egglib
        multiqc primer3
        python2 wasp test'''.split()
    envs2 = '''base snk work'''.split()
    envs = envs1
    if args.opt == 2:
        envs = envs2
    print("will update %d environments: %s" % (len(envs), ' '.join(envs)))
    for env in envs:
        print("updating %s" % env)
        p = Popen(["conda update -n %s --all" % env], stdin=PIPE, shell=True)
        outs, errs = p.communicate(input=b'y\n')
        p.terminate()
        sh("conda env export -n %s --no-builds > $snk/envs/%s.yml" % (env, env))

def merge_dirs(args):
    diris, diro = args.diri, args.diro
    mkdir(diro, overwrite=True)
    for diri in args.diri:
        for fn in os.listdir(diri):
            fi = op.join(diri, fn)
            fo = op.join(diro, fn)
            if not op.isfile(fi): continue
            if op.isfile(fo):
                if not cmp(fi, fo):
                    if args.replace:
                        copy(fi, fo)
                    else:
                        print("%s/%s diff from %s/%s - skipped" % (diri, fn, diro, fn))
            else:
                copy(fi, fo)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'common bash commands'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("update_conda",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "update selected conda environments")
    sp1.add_argument('--opt', type=int, default=1, help = 'option')
    sp1.set_defaults(func = update_conda)

    sp1 = sp.add_parser("merge_dirs",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "merge two or more directories into a single one")
    sp1.add_argument('diro', help = 'output directory')
    sp1.add_argument('diri', nargs='+', help = 'one or more input directories')
    sp1.add_argument('--replace', action='store_true', help = 'replace file?')
    sp1.set_defaults(func = merge_dirs)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


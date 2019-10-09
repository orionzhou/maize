#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from subprocess import Popen, PIPE, run
from maize.apps.base import sh

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

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


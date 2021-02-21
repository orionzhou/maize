#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from subprocess import Popen, PIPE, run
from jcvi.apps.base import sh, mkdir

def pull(args):
    repos = args.repos.split(" ")
    for repo in repos:
        dir1 = op.join(os.getenv("HOME"), 'git', repo)
        dir2 = op.join(os.getenv("HOME"), 'projects', repo)
        if op.isdir(dir1):
            os.chdir(dir1)
            logging.debug(f"repo {repo}: {dir1}")
        elif op.isdir(dir2):
            os.chdir(dir2)
            logging.debug(f"repo {repo}: {dir2}")
        else:
            logging.error(f"repo {repo}: not found - skipped")
        sh("git stash", log=False)
        sh("git pull", log=False)
        sh("git stash pop", log=False)
        sh("git commit -am 'normal commit'", log=False)
        sh("git push origin master", log=False)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'common bash commands'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    repos = "appconfig maize nf rmaize demo assets atlas barn biomap bsseq cage chipseq cold cre epi genome grn misc ml reseq rnaseq s3 stress wgc"
    repos2 = "maizeumn.github.io orionzhou.github.io nf-core-methylseq nf-core-rnaseq nf-core-sarek nf-core-chipseq"
    sp1 = sp.add_parser("pull",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "pull/push updates from/to multiple git repos")
    sp1.add_argument('--repos', default=repos, help = 'git repos')
    sp1.set_defaults(func = pull)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


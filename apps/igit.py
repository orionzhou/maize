#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

import git
from git import Repo
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

def main(args):
    repos = args.repos.split(" ")
    commit_msg = 'normal commit'
    for repo in repos:
        dir1 = op.join(os.getenv("HOME"), 'git', repo)
        dir2 = op.join(os.getenv("HOME"), 'projects', repo)
        diri = dir1
        if op.isdir(dir1):
            logging.info(f"repo {repo}: {dir1}")
        elif op.isdir(dir2):
            diri = dir2
            logging.info(f"repo {repo}: {dir2}")
        else:
            logging.error(f"repo {repo}: not found - skipped")
            continue
        repo = Repo(diri)
        push = repo.is_dirty()
        if push:
            repo.git.stash("save")
        repo.remotes.origin.pull()
        if push:
            repo.git.stash("pop")
            repo.git.add(update=True)
            repo.index.commit(commit_msg)
            origin = repo.remote(name='origin')
            origin.push()

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
         formatter_class = argparse.ArgumentDefaultsHelpFormatter,
         description = 'sync multiple git repos'
    )
    repos = "appconfig maize nf rmaize demo assets atlas barn biomap bsseq cage chipseq cre epi genome grn misc ml reseq rnaseq s3 stress wgc"
    repos2 = "maizeumn.github.io orionzhou.github.io nf-core-methylseq nf-core-rnaseq nf-core-sarek nf-core-chipseq"
    repos = 'maize cre'
    ps.add_argument('--repos', default=repos, help = 'git repos')
    args = ps.parse_args()
    main(args)

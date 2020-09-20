#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
#import boto3

from maize.apps.base import sh, mkdir

def s3_sync(args):
    os.chdir(args.dirh)
    bucket = "%s%s" % (args.pre, args.bucket)
    dry = "--dry-run" if args.dry else ''
    cmd = "s3cmd -c %s/appconfig/s3cfg2" % os.environ['git'] if args.personal else "s3cmd"
    if args.dry:
        sh("%s sync %s --delete-removed %s/ s3://%s/" % (cmd, dry, bucket, bucket))
    else:
        sh("%s sync %s --delete-removed %s/ s3://%s/" % (cmd, dry, bucket, bucket))
        sh("%s setacl -P -r s3://%s" % (cmd, bucket))


if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Amazon S3 utilities'
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("sync",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "sync an s3 bucket")
    sp1.add_argument('bucket', choices=['igv','igv-data','multiqc'], help = 'bucket suffix')
    sp1.add_argument('--dirh', default="%s/s3" % os.environ["proj"], help = 'local directory for S3')
    sp1.add_argument('--pre', default="zhoup-", help="bucket prefix")
    sp1.add_argument('--personal', action="store_true", help="use personal aws account instead of msi account?")
    sp1.add_argument('--dry', action="store_true", help="dry run?")
    sp1.set_defaults(func = s3_sync)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()


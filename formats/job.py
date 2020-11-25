#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from jinja2 import Template
from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

tmpl_slurm = '''#!/bin/bash -l
#SBATCH -N {{ nodes }} -n {{ ntasks }} -c {{ cores }}
#SBATCH --time={{ time }}
#SBATCH --mem={{ mem }}
#SBATCH -p {{ partitions }}
#SBATCH --mail-type={{ mail_type }}
#SBATCH --mail-user={{ email }}
#SBATCH -o {{ out }}

echo ${SLURM_JOBID}
'''

def create_job_script(args):
    tmpl = ''
    if args.fmt == 'slurm':
        tmpl = Template(tmpl_slurm)
    else:
        print("unsupported format: %s" % args.fmt)
        sys.exit(1)
    msg = tmpl.render(nodes = args.nodes,
                      ntasks = args.ntasks,
                      cores = args.cores,
                      time = args.time,
                      mem = args.mem,
                      partitions = args.partitions,
                      mail_type = args.mail_type,
                      email = args.email,
                      out = args.job_out
    )

    fho = must_open(args.out, 'w')
    fho.write(msg)
    fho.close()

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "create (slurm) job script from templates"
    )

    ps.add_argument('--out', '-o', default='job', help = 'output job script')
    ps.add_argument('--fmt', default='slurm', help = 'job scheduler')
    ps.add_argument('--nodes', '-N', default=1, help='num. nodes')
    ps.add_argument('--ntasks', '-n', default=1, help='num. tasks')
    ps.add_argument('--cores', '-c', default=1, help='num. cores')
    ps.add_argument('--time', '-t', default="10:00:00", help='wall time')
    ps.add_argument('--mem', '-m', default="20gb", help='memory')
    ps.add_argument('--mail_type', default="FAIL", help='when to send email')
    ps.add_argument('--email', default="zhoux379@umn.edu", help='email')
    ps.add_argument('--partitions', '-p', default="small,amdsmall", help='partition/queue')
    ps.add_argument('--job_out', default="%x.out", help='job output')

    args = ps.parse_args()
    create_job_script(args)


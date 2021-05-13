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
#SBATCH -p {{ queue }}
#SBATCH --mail-type={{ mail_type }}
#SBATCH --mail-user={{ email }}
#SBATCH -o {{ out }}

echo ${SLURM_JOBID}
'''

tmpl_sge = '''#!/bin/bash -l
#$ -N {{ name }}
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -q {{ queue }}
#$ -l h_vmem={{ mem }}
#$ -l s_rt={{ time }}
#$ -pe smp {{ ntasks }}
#$ -l h_vmem={{ mem }}
#$ -m {{ mail_type }}
#$ -M {{ email }}
#$ -j y
#$ -o {{ out }}
'''


def create_job_script(args):
    tmpl = ''
    nodes = args.nodes
    ntasks = args.ntasks
    cores = args.cores
    time = args.time
    mem = args.mem
    queue = args.queue
    mail_type = args.mail_type
    email = args.email
    name = args.job_name
    out = args.job_out
    if args.fmt == 'slurm':
        tmpl = Template(tmpl_slurm)
        queue = 'small,amdsmall'
        mail_type = "FAIL"
        mem = '20gb'
    elif args.fmt == 'sge':
        tmpl = Template(tmpl_sge)
    else:
        print("unsupported format: %s" % args.fmt)
        sys.exit(1)
    msg = tmpl.render(nodes = nodes,
        ntasks = ntasks,
        cores = cores,
        time = time,
        mem = mem,
        queue = queue,
        mail_type = mail_type,
        email = email,
        name = name,
        out = out
    )

    fho = must_open(args.out, 'w')
    fho.write(msg)
    fho.close()

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "create job script (slurm/sge) from templates"
    )

    fmts = ['slurm','sge']
    ps.add_argument('--out', '-o', default='job', help = 'output job script')
    ps.add_argument('--fmt', default='sge', choices=fmts, help = 'job scheduler')
    #ps.add_argument('--queue', '-q', default="smallfat.q", help='queue/partition')
    ps.add_argument('--queue', '-q', default="all.q", help='queue/partition')
    ps.add_argument('--nodes', '-N', default=1, help='num. nodes')
    ps.add_argument('--ntasks', '-n', default=1, help='num. tasks')
    ps.add_argument('--cores', '-c', default=1, help='num. cores')
    ps.add_argument('--time', '-t', default="10:00:00", help='wall time')
    ps.add_argument('--mem', '-m', default="20G", help='memory')
    ps.add_argument('--mail_type', default="a", help='when to send email')
    ps.add_argument('--email', default="zhoupeng03@caas.cn", help='email')
    ps.add_argument('--job_name', default="test", help='job name')
    ps.add_argument('--job_out', default="%x.out", help='job output')

    args = ps.parse_args()
    create_job_script(args)


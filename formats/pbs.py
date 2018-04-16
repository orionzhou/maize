#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from string import Template

from maize.formats.base import must_open

class subPbsJob(object):
    
    def __init__(self, cmds = [], prefix = 'subjob', \
            opt = 'bash', thread = 1, \
            bash = 'bash', parallel = 'parallel'):
        assert isinstance(cmds, list), "cmds must be a list: %s" % str(cmds)
        self.fname = "%s.sh" % prefix
        self.cmds = cmds
        self.bash = bash
        self.parallel = parallel
        self.jcmd = ''
        if opt == 'bash':
            self.jcmd = "%s %s" % (self.bash, self.fname)
        elif opt == 'parallel':
            self.jcmd = "%s -j %d < %s" % (self.parallel, int(thread), self.fname)
        else:
            logging.error("unknown option: %s" % opt)
            sys.exit(1)

    def add_cmd(self, cmd):
        if isinstance(cmd, str):
            self.cmds.append(cmd)
        elif isinstance(cmd, list):
            self.cmds += cmd
        else:
            logging.error("unknown cmd: ", cmd)

    def write(self):
        fho = must_open(self.fname, "w")
        fho.write("\n".join([x.strip() for x in self.cmds]) + "\n")
        fho.close()

    def __str__(self):
        return "Bash script [%s, %s] with %d cmd-lines" % \
            (self.fname, self.opt, len(self.cmds))

    __repr__ = __str__

class PbsJob(object):
    __stanza__ = '''#PBS -q $queue
#PBS -l nodes=$node:ppn=$ppn
#PBS -l walltime=$walltime
#PBS -m ae
#PBS -M $email

$cmds'''
    
    def __init__(self, queue = 'small', node = 1, ppn = 1, 
            walltime = "10:00:00", 
            email = 'zhoux379@umn.edu', 
            prefix = "jobpre",
            parallel = "parallel",
            bash = "bash"):
        self.queue = queue
        self.node = node
        self.ppn = ppn
        self.walltime = walltime
        self.email = email
        self.prefix = prefix
        self.parallel = parallel
        self.bash = bash
        self.fname = "%s.pbs" % prefix
        self.subjobs = []
        self.cmds = []
    
    @classmethod
    def from_cfg(cls, jcfg = {}, jcmds = [], bcfgs = [], bcmds = [],
            njob = 0, prefix = 'jobpre', bash = 'bash', parallel = 'parallel'):
        queue = jcfg.get("queue", "small")
        node = jcfg.get("node", 1)
        ppn = jcfg.get('ppn', 1)
        walltime = jcfg.get("walltime", "10:00:00")
        email = jcfg.get("email", "zhoux379@umn.edu")
        job = cls(queue = queue,
                node = node,
                ppn = ppn,
                walltime = walltime,
                email = email,
                prefix = prefix,
                bash = bash,
                parallel = parallel
        )
        job.add_cmd(jcmds)

        if njob > 0:
            if bcfgs is None or len(bcfgs) == 0:
                bcfgs = [{} for i in range(njob)]
            if bcmds is None or len(bcmds) == 0:
                bcmds = [[] for i in range(njob)]
        assert njob == len(bcmds) == len(bcfgs), "not %d subjobs: %s" % (njob, str(bcmds))
        for j in range(njob):
            bcfg = bcfgs[j]
            opt = bcfg.get("opt", "bash")
            thread = bcfg.get("thread", 1)
            job.add_subjob(cmds = bcmds[j], opt = opt, thread = thread)
        return job
 
    def add_subjob(self, cmds = [], opt = 'bash', thread = None):
        jpre = "%s.%d" % (self.prefix, len(self.subjobs) + 1)
        if thread is None: thread = self.ppn
        subjob = subPbsJob(cmds = cmds, prefix = jpre, \
                opt = opt, thread = thread, \
                bash = self.bash, parallel = self.parallel)
        self.subjobs.append(subjob)
        self.cmds.append(subjob.jcmd)
    
    def add_cmd(self, cmd):
        if isinstance(cmd, str):
            self.cmds.append(cmd)
        elif isinstance(cmd, list):
            self.cmds += cmd
        else:
            logging.error("unknown cmd: ", cmd)

    def write(self):
        for subjob in self.subjobs:
            subjob.write()
        source = Template(self.__stanza__)
        cmd_str = "\n".join([x.strip() for x in self.cmds])
        params = {
                "queue": self.queue,
                "node": self.node,
                "ppn": self.ppn,
                "walltime": self.walltime,
                "email": self.email,
                "cmds": cmd_str + "\n"
        }
        fhj = must_open(self.fname, "w")
        fhj.write(source.substitute(params))
        fhj.close()

    def __str__(self):
        return "PBS script [%s] with %d sub-jobs" % (self.fname, len(self.subjobs))

    __repr__ = __str__

def create_job_chain(fjs, fo):
    cmds = ['#!/bin/bash'] 
    jobs = []
    for i in range(len(fjs)):
        fj = fjs[i]
        assert op.isfile(fj), "cannot read %s" % fj
        job = "job%d" % (i+1)
        if i == 0:
            cmds.append("%s=$(qsub %s)" % (job, fj))
            cmds.append("echo $%s" % job)
        else:
            pjob = jobs[i-1]
            cmds.append("%s=$(qsub -W depend=afterok:$%s %s)" % (job, pjob, fj))
            cmds.append("echo $%s" % job)
        jobs.append(job)
    fho = must_open(fo, "w")
    fho.write("\n".join(cmds)+"\n")

def start_pbs_jobs(fjs):
    import subprocess
    import re
    ptnj = re.compile("^(\d+)\.([\w\.]+)$")
    jids = []
    for i in range(len(fjs)):
        fj = fjs[i]
        assert op.isfile(fj), "cannot read %s" % fj
        if i == 0:
            proc = subprocess.Popen(["qsub", fj], stdout=subprocess.PIPE)
            out = str(proc.communicate()[0], 'utf-8').strip()
        else:
            pjid = jids[i-1]
            argu = "-W depend=afterok:%s" % pjid
            proc = subprocess.Popen(["qsub", argu, fj], stdout=subprocess.PIPE)
            out = str(proc.communicate()[0], 'utf-8').strip()
        res = out.match(ptnj)
        if res:
            jid = res.group(1)
            jids.append(jid)
        else:
            eprint("qsub error: %s" % out)
            sys.exit(1)
    print("\n".join(jids))

if __name__ == '__main__':
    x = PbsJob()
    x.add_subjob(['ls -l'], opt = 'parallel')
    print(x, x.queue)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
from string import Template

from maize.formats.base import must_open

class PbsJob(object):
    __stanza__ = \
"""#PBS -l nodes=1:ppn=$ppn,walltime=$walltime
#PBS -m ae
#PBS -M $email
#PBS -q $queue

$cmds
"""
    
    def __init__(self, ppn = 1, 
            walltime = "2:00:00", 
            email = 'zhoux379@umn.edu', 
            queue = 'small',
            cmds = "echo $HOSTNAME"):
        self.params = {
                "queue": queue,
                "walltime": walltime,
                "ppn": ppn,
                "email": email,
                "cmds": cmds
        }
        self.source = Template(self.__stanza__)

    def __str__(self):
        return self.source.substitute(self.params)

    __repr__ = __str__

    def write(self, outfile):
        fh = must_open(outfile, "w")
        fh.write(self.source.substitute(self.params))
        fh.close()

if __name__ == '__main__':
    print(PbsJob())

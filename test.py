# -*- coding: utf-8 -*-
import sys
import os.path
sys.path.append(os.environ['PYTHON_HOME'] + '/lib/python2.6/site-packages')
sys.path.append(os.environ['PYTHON_HOME'] + '/lib64/python2.6/site-packages')
sys.path.append(os.environ['SCRIPT_HOME_PYTHON'])

from time import clock, time
from InitPath import *
from Seq import *

if __name__ == "__main__":
    seq_db = "/project/youngn/zhoup/Data/misc3/spada.crp/Athaliana/01_genome/01_refseq.fa"
    seq_id = "Chr2"
    itv = 100000

    for i in range(1, 10):
        seq_beg = itv * i + 1
        seq_end = itv * i + 10
        
        c0, t0 = clock(), time()
        print seqret(seq_db, seq_id, seq_beg, seq_end, -1)
        
        print "process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0)
        

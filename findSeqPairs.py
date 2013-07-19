# -*- coding: utf-8 -*-
import sys
import os.path
sys.path.append(os.environ['PYTHON_HOME'] + '/lib/python2.6/site-packages')
sys.path.append(os.environ['PYTHON_HOME'] + '/lib64/python2.6/site-packages')
sys.path.append(os.environ['SCRIPT_HOME_PYTHON'])

from InitPath import *
from Common import *
from Seq import *
import tablet as T

if __name__ == "__main__":
    dir = os.path.join(DIR_Misc3, "spada.crp/Athaliana")
    f_seq = os.path.join(dir, "01_genome/01_refseq.fa")
    fi = os.path.join(dir, "31_model_SPADA/61_final.tbl")
    fo = os.path.join(dir, "seq_nt.tbl")
    
    t1 = T.read(fi, delim="\t")
    fho = open(fo, "w")

    for row in t1:
        id, fam, chr, beg, end, str, e, seq = getVar(row, range(7) + [10])
    
        seq_str = seqret(f_seq, chr, int(beg), int(end), str)

        fho.write("-".join([id, fam, chr, beg, end, seq_str]))

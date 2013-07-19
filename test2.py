import sys
import os.path
sys.path.append(os.environ['PYTHON_HOME'] + '/lib/python2.6/site-packages')
sys.path.append(os.environ['PYTHON_HOME'] + '/lib64/python2.6/site-packages')
sys.path.append(os.environ['SCRIPT_HOME_PYTHON'])
from time import clock, time
from subprocess import Popen, PIPE
from Bio import SeqIO

def seqret(f_seq, id, beg, end, strand=1):
    if not os.path.isfile(f_seq):
        sys.exit(f_seq + " is not there")
    f_idx = f_seq + ".idx"
    seq_idx = SeqIO.index_db(f_idx, f_seq, "fasta")
    if not id in seq_idx:
        sys.exit("cannot find " + id + "in " + f_seq)
    
    seq_obj = seq_idx[id].seq[beg-1:end]
    if strand == -1 or strand == "-":
        seq_obj = seq_obj.reverse_complement()
    seq_str = seq_obj.tostring()
    return seq_str

if __name__ == "__main__":
   
    seq_db = "/project/youngn/zhoup/Data/misc3/spada.crp/Athaliana/01_genome/01_refseq.fa"
    seq_id = "Chr3"
    itv = 100000
    seq_beg = itv * 1 + 1
    seq_end = itv * 1 + 10
        
    c0, t0 = clock(), time()
    print seqret(seq_db, seq_id, seq_beg, seq_end, 1)

    print "process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0)

    c0, t0 = clock(), time()
    f_perl = "/soft/bioperl/5.16.1/bin/bp_seqret.pl"
    cmd = '%s -d %s -i %s -s %d -e %d' % (f_perl, seq_db, seq_id, seq_beg, seq_end)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    lines = p.stdout.readlines()
    print lines
    print "process time [%.2fs], wall time [%.2fs]" % (time()-t0, clock()-c0)
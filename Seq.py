# -*- coding: utf-8 -*-
import sys
import os.path
sys.path.append(os.environ['PYTHON_HOME'] + '/lib/python2.6/site-packages')
sys.path.append(os.environ['PYTHON_HOME'] + '/lib64/python2.6/site-packages')
sys.path.append(os.environ['SCRIPT_HOME_PYTHON'])
from Bio import Seq, SeqIO
from subprocess import Popen, PIPE

def rmSeqDesc(fi, fo):
    seqRecdLst = [];
    with open(fi, "r") as fhi:
        for seqRecd in SeqIO.parse(fhi, "fasta"):
            print "%s\t%s" % (seqRecd.id, seqRecd.description);
            seqRecd.description = "";
            seqRecdLst.append(seqRecd);
    with open(fo, "w") as fho:
        SeqIO.write(seqRecdLst, fho, "fasta");

def seqret(f_seq, seq_id, seq_beg, seq_end, strand=1):
    f_perl = "/soft/bioperl/5.16.1/bin/bp_seqret.pl"
    cmd = '%s -d %s -i %s -s %d -e %d' % (f_perl, f_seq, seq_id, seq_beg, seq_end)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    
    for line in p.stderr:
        print "!!!!! error retrieving sequence !!!!!\n" + line
        sys.exit(1)

    id, seqstr = '', ''
    for line in p.stdout:
        line = line.rstrip()
        if line.startswith('>'):
            id = line
        elif line.isalpha():
            seqstr += line
    
    if strand == -1 or strand == "-":
        seqstr = Seq.Seq(seqstr).reverse_complement().tostring()
    return seqstr

def seqret_slow(f_seq, id, beg, end, strand=1):
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
    fi = "/project/youngn/zhoup/Data/misc3/spada.crp/Athaliana/01_genome/01_refseq.fa"
    print seqret(fi, 'Chr1', 100001, 100010, -1)

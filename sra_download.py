#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: sra_download
    :platform: Unix, Windows
    :synopsis: A python script to download SRA files and extract Fastq from
        them
    
.. moduleauthor:: Peng Zhou <zhoupenggeni@gmail.com>

"""

import os
import os.path as op
import numpy as np
from string import Template
import sys
import time
import argparse

def main():
    """This function downloads SRA files from NCBI and extract Fastq files
    
    Args:
        dirw (str): The working directory (will be created if not already
            exists)

        sralist (str): A tabular seperated file (\*.tsv) specifying the SRA
            IDs ("SampleID", required) and meta infomation (optional).
            Sample sralist file:

==========   ========   ===============  ========  ========= 
SampleID     Species    Tissue           Genotype  Treatment
==========   ========   ===============  ========  ========= 
SRR957415    Zea mays   Endosperm_12DAP  B73
SRR957420    Zea mays   Pericarp         B73
SRR957421    Zea mays   Endosperm_Crown  B73
==========   ========   ===============  ========  =========

    Kwargs:
        qscript (str): PBS job template used to create the job script
        
    Output:
        * Sub-folder `03.sra` contains SRA files (\*.sra) downloaded
            from NCBI
        * Job script `04.fastqdump.pbs` that can be submitted to MSI queue 
            (e.g., mesabi) for executation
            
            - Make manual modifications if necessary (at least remember to 
                change Email address where notifications will be sent
            - Submit job by ``qsub 04.fastqdump.pbs`` 

        * Sub-folder `05.reads` contains Fastq files (\*.1.fastq,
            \*.2.fastq)

    Raises:
        AttributeError, KeyError

    To see help information:

.. code-block:: console

    usage: sra_download.py [-h] dirw sralist qscript

    Download SRA files and extract fastq sequences

    positional arguments: 
      dirw        working directory (default: 
                  /scratch.global/zhoux379/shortread/test) 
      sralist     list file of SRA IDs (default: 00.0.srr.tsv)

    optional arguments: 
      -h, --help  show this help message and exit
      qscript     PBS script template (default:
                  /home/springer/zhoux379/git/robin/templates/sra_download.pbs)

    """
    parser = argparse.ArgumentParser(
        description = 'Download SRA files and extract fastq sequences'
    )
    parser.add_argument(
        'dirw', default = "/scratch.global/zhoux379/shortread/test", help = 'working directory (default: /scratch.global/zhoux379/shortread/test)'
    )
    parser.add_argument(
        'sralist', default = "00.0.srr.tsv", help = 'list file of SRA IDs (default: 00.0.srr.tsv)'
    )
    parser.add_argument(
        '--qscript', default = "/home/springer/zhoux379/git/robin/templates/sra_download.pbs", help = 'PBS script template (default: /home/springer/zhoux379/git/robin/templates/sra_download.pbs)'
    )
    args = parser.parse_args()

    dirw, sralist, ftem = args.dirw, args.sralist, args.qscript
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(sralist), "%s is not a file" % sralist 
    ary = np.genfromtxt(sralist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    if not op.isdir("03.sra"): os.makedirs("03.sra")
    if not op.isdir("05.reads"): os.makedirs("05.reads")
    d03 = "03.sra" #op.abspath("03.sra")
    f04 = "04.fastqdump.sh" #op.abspath("04.fastqdump.sh")
    d05 = "05.reads" #op.abspath("05.reads")
    fhc = open(f04, "w")
    os.chdir("03.sra")
    fnames = []
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        rid = row[0]
        print(rid)
        cmd = "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra" % (rid[0:3], rid[0:6], rid, rid)
        #os.system(cmd)
        f_sra = "%s.sra" % rid
        assert op.isfile(f_sra), "%s is not there" % f_sra
        fnames.append(f_sra)
        fhc.write("fastq-dump --gzip --split-files -outdir %s %s/%s\n" % (d05, d03, f_sra))
    fname_str = " ".join(fnames)
    fhc.close()
    os.chdir("..")

    assert op.isfile(ftem), "cannot read template file:  %s" % ftem
    fht = open(ftem, "r")
    src = Template(fht.read())

    fjob = "04.fastqdump.pbs"
    fhj = open(fjob, "w")
    fhj.write(src.substitute({'DIR':op.abspath(dirw)}))
    fhj.close()
    #os.system("qsub %s -v DIR=%s" % (name, dirw))

if __name__ == "__main__":
    main()

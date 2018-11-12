#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import argparse

def read_cfg(fc):
    (orgs, fis) = ([], [])
    fhc = open(fc, "r")
    for line in fhc:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        (org, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            print os.access(fi, os.F_OK)
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
    fhc.close()
    return (orgs, fis)
def orthomcl_adjustfasta(orgs, fis, diro):
    if not op.exists(diro): os.makedirs(diro)
    os.chdir(diro)
    for i in range(0,len(orgs)):
        (org, fi) = (orgs[i], fis[i])
        cmd = "orthomclAdjustFasta %s %s 1" % (org, fi)
        print cmd
        os.system(cmd)
    os.chdir("..")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        'cfgfile', help = 'config file (.csv)'
    )
    parser.add_argument(
        'outdir', help = 'output directory'
    )
    args = parser.parse_args()
    (fc, dirw) = (args.cfgfile, args.outdir)
    dirw = op.realpath(dirw)
    (orgs, fis) = read_cfg(fc)
    
    if not op.exists(dirw): os.makedirs(dirw)
    
    cdir = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] = os.environ['PATH']+':'+cdir
    cwd = os.getcwd()
    os.chdir(dirw)
    
    #orthomcl_adjustfasta(orgs, fis, "01_fasta")
    #os.system("orthomclFilterFasta 01_fasta 10 20")
    #os.system("makeblastdb -dbtype prot -in goodProteins.fasta -out allprot")
    #os.system("qsub.blat.pl -i goodProteins.fasta -o 11.blast -t %s/allprot -p 2" % dirw)
    
    #os.system("cat 11.blast/*.tbl > 12.blast.tbl")
    #os.system("orthomclBlastParser 12.blast.tbl 01_fasta > 15.sim.seq.txt")
    #os.system("orthomclInstallSchema orthomcl.config log.install")
    #os.system("orthomclLoadBlast orthomcl.config 15.sim.seq.txt")
    os.system("orthomclPairs orthomcl.config log.run cleanup=no startAfter=useLog")
    #os.system("orthomclDumpPairsFiles orthomcl.config")
    #os.system("mcl mclInput --abs -I 1.5 -o 21.mcl.out")
    #os.system("orthomclMclToGroups mt 1000 < 21.mcl.out > 25.groups.txt")
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
from igv import IGV
igv = IGV()

DIR_data = "/home/youngn/zhoup/Data"
#dir = "%s/%s/%s" % (DIR_data, "misc4", "spada.crp.Athaliana")
dir = "%s/%s/%s" % (DIR_data, "misc4", "spada.defl.Mtruncatula_3.5")

fi = "%s/%s" % (dir, "31_model_evaluation/91_compare.tbl")
dirO = "%s/%s" % (dir, "42_imgs")

#org = "Athaliana"
org = "Mtruncatula_3.5"
igv.genome(org)

f_gen = "%s/misc4/genome/%s/51_gene.gff" % (DIR_data, org)
#f_gs = "%s/misc2/crp.gs/crp.at.gff" % (DIR_data)
f_gs = "%s/misc2/crp.gs/crp.mt35.gff" % (DIR_data)
f_spa_all = "%s/21_model_prediction/26_all.gff" % (dir)
f_spa = "%s/31_model_evaluation/61_final.gff" % (dir)
f_bam = "%s/misc2/spada.rnaseq/%s/21_tophat/accepted_hits.bam" % (DIR_data, org)

igv.load(f_gen)
igv.load(f_gs)
igv.load(f_spa_all)
igv.load(f_spa)
igv.load(f_bam)
igv.send("squish")

if __name__ == "__main__":
    fhi = open(fi, "r");
    if not os.path.exists(dirO):
        os.makedirs(dirO)
    for line in fhi:
        line = line.strip("\n")
        ps = line.split("\t")
        if line == "" or ps[0] == "idQ":
            continue;
        idQ = ps[0]
        tag = int(ps[2])
        fo = "%s/%s.png" % (dirO, idQ)
        if(tag > 2):
            print idQ, tag
            igv.go(idQ)
            igv.save(fo)

	

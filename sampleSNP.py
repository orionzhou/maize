# -*- coding: utf-8 -*-
import my;
import sys;
import re;
import os.path;
import random;

fIn = os.path.join(my.DIR_data, "Variant_Mt");
fInfo = os.path.join(my.DIR_data, "Variant_Mt", "variants.info");
lstAcc = ["HM001","HM002","HM003","HM004","HM005_gsnapout8mm", \
    "HM006","HM009","HM011","HM015","HM017-I", \
    "HM029_gsnapout8mm", "HM101_gsnapout8mm"];
lstAcc = ["HM001","HM002","HM003","HM004","HM005_gsnapout8mm", \
    "HM006","HM009","HM011","HM015", \
     "HM101_gsnapout8mm"];

numberRandom = 2000;

if __name__ == "__main__":   
    randDict = dict();
    with open(fInfo, "r") as fInfoH:
        for line in fInfoH:
            rst = re.match("^MtChr(\d)\t(\d+).*", line);
            if rst != None:
                chr = int(rst.group(1));
                cntLine = int(rst.group(2));
                randLst = random.sample(range(1,cntLine), numberRandom);
                randDict[chr] = randLst;
    myRet = my.VariantRet(lstAcc);
    
    randDict = {2:randDict[2]};
    for chr,randLst in sorted(randDict.iteritems()):
        fOut = os.path.join(my.DIR_output, "SNP_MtChr%d.txt"%chr);
        myRet.retByLine(chr, randLst, fOut);

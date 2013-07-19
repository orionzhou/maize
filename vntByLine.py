# -*- coding: utf-8 -*-
import my;
import sys;
import re;
import os.path;
import random;

lstAcc = ["HM001","HM002","HM003","HM004","HM005_gsnapout8mm", \
    "HM006","HM009","HM011","HM015","HM017-I", \
    "HM029_gsnapout8mm", "HM101_gsnapout8mm"];
chrNameLst = ["MtChloro", "MtMito"] + ["MtChr%d" % i for i in range(1,9)];
chrNameLst = [chrNameLst[0]];

fVnt = os.path.join(my.DIR_Out, "vnt_{0}.txt");
fReport = os.path.join(my.DIR_Stat, "stat_{0}.txt");

def getRandomLine(chrName, numberRandom=10, vntType="snp"):
    fInfo = os.path.join(my.DIR_data, "Variant_Mt", "variants.info");
    flag = False;
    with open(fInfo, "r") as fInfoH:
        for line in fInfoH:
            line = line.strip("\n");
            if line == "":
                break;
            ele = line.split("\t");
            chrNameInFile = ele[0];
            if chrName == chrNameInFile:
                assert vntType in ("snp", "indel");
                if vntType == "snp":
                    cntLine = int(ele[2]);
                else:
                    cntLine = int(ele[3]);
                randLst = random.sample(range(1,cntLine), numberRandom);
                flag = True;
                print "%d lines of %s(s) out of %d sampled" % (numberRandom, vntType, cntLine);
    assert flag == True;

if __name__ == "__main__":   
    myRet = my.VntRet(lstAcc);
    number = 10;
    for chrName in chrNameLst:
        #lineLst = getRandomLine(chrName,numberRandom);
        lineLst = range(1,10000);
        myRet.retByLine(chrName, lineLst, fVnt.format("%s_%d"%(chrName,number));
        myVntOut = my.VntReport(cutOffCov, cutOffFreq, cutOffUniq);
        myVntOut.readIn(fVnt, fReport.format("%s_%d"%(chrName,number));
    

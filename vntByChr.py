# -*- coding: utf-8 -*-
import os.path;
import my;

lstAcc = ["HM001","HM002","HM003","HM004","HM005_gsnapout8mm", \
    "HM006","HM009","HM011","HM015","HM017-I", \
    'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];

cutOffCov, cutOffFreq, cutOffUniq= 2, 0.7, 2;
suffix1 = "chr10k";
suffix2 = "snp";
#suffix2 = "indel";
fIn = os.path.join(my.DIR_In, "location_{0}").format(suffix1);
assert os.path.exists(fIn);
DIR_Work = os.path.join(my.DIR_Misc, "Vnt_"+suffix);
if not os.path.exists(DIR_Work):
    os.mkdir(DIR_Work);
fVnt = os.path.join(DIR_Work, "vnt_{0}_{1}.txt");
fReport = os.path.join(DIR_Work, "vntStat_{0}_{1}.txt");

if __name__ == "__main__":   
    locDict = my.getLoc(fIn, 2);
    for chrName in locDict:
        chrDict = {chrName:locDict[chrName]};
        fVnt = fVnt.format(suffix2, chrName);
        fReport = fReport.format(suffix2, chrName);
        print "\n".join([fVnt, fReport]);
        myVntRet = my.VntRet(lstAcc);
        myVntRet.retByLoc(chrDict, fVnt, suffix2);
        myVntOut = my.VntReport(cutOffCov, cutOffFreq, cutOffUniq);
        myVntOut.readIn(fVnt, fReport);

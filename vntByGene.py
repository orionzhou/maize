# -*- coding: utf-8 -*-
import my;
import os.path;

lstAcc = ["HM001","HM002","HM003","HM004","HM005_gsnapout8mm", \
    "HM006","HM009","HM011","HM015","HM017-I", \
    'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];
lstAcc = ["HM005_gsnapout8mm", 'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];

cutOffCov, cutOffFreq, cutOffUniq= 2, 0.7, 2;

suffix1 = "DEFL";
suffix2 = "snp";
#suffix2 = "indel";
fIn = os.path.join(my.DIR_In, "location_{0}").format(suffix1);
assert os.path.exists(fIn);
DIR_Work = os.path.join(my.DIR_Misc, suffix1);
if not os.path.exists(DIR_Work):
    os.mkdir(DIR_Work);
fVnt = os.path.join(DIR_Work, "vnt_{0}.txt".format(suffix2));
fReport = os.path.join(DIR_Work, "vntStat_{0}.txt".format(suffix2));
#fVntDesc = os.path.join(DIR_Work, "vntDesc_{0}.txt".format(suffix2));

if __name__ == "__main__":   
    locDict = my.getLoc(fIn, 2, "F", 3); #only read 3+ columns
    myVntRet = my.VntRet(lstAcc);
    myVntRet.retByLoc(locDict, fVnt, suffix2);
    #myVntOut = my.VntReport(cutOffCov, cutOffFreq, cutOffUniq, "gene", fIn);
    #myVntOut.readIn(fVnt, fReport);

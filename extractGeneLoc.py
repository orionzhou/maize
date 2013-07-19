# -*- coding: utf-8 -*-
import sys, os.path, re;
import my;

fIn = os.path.join(my.DIR_input, "blastLoc.txt");
fOut = os.path.join(my.DIR_input, "geneLoc_CCP.txt");

geneFamily = "CCP";
if __name__ == "__main__":
    #get gene location for a gene family
    fInH = open(fIn, "r");
    fOutH = open(fOut, "w");
    cntGene = 0;
    posStrDict = dict();
    for line in fInH:
        line = line.strip("\n");
        if line == "":
            break;
        rst = re.findall("([A-Z]{2,3}[\d\.]+)", line);
        if rst is not None:
            assert len(rst) == 1;
            geneId = rst[0];
            rst = re.findall("(MtChr\d\:\d+\.\.\d+)", line);
            assert len(rst) >= 1;
            posStrDict[geneId] = ";".join(rst);
    print >>fOutH, ">%s(%d)" % (geneFamily, len(posStrDict.keys()));
    for geneId, posStr in sorted(posStrDict.iteritems()):
        print >>fOutH, "\t".join([geneId, "", posStr, ""]);
        

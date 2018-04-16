# -*- coding: utf-8 -*-
import sys
import os
import re

def getLoc(fIn, option1=1, option2="F", option3=1, fOut=None):
# generate location dict 
# option1 = 1/2 hierarchical/chromosomal
# option2 = 'F'/'R' forward/reverse dictionary
# option3 = 1/2/3 to get all/CDS+Intron+UTR/Gene columns 
    fInH = open(fIn, "r");
    type = fInH.readline().strip("\n");
    assert type in ["Gene", "Chromosome", "Probe"];
    if option3 == 2 or option3 == 3:
        assert type == "Gene";
    if fOut != None:
        fOutH = open(fOut, "w");
    locDict = dict();
    if option1 == 2:
        for line in fInH:
            line = line.strip("\n");
            if line == "":
                continue;
            ary = line.split("\t",2);
            if option3 == 2:
                line = ary[-1];
            elif option3 == 3:
                line = ary[1];
            pattern = re.compile(\
                "(MtChr[1-8]|MtChloro|MtMito)[:\-\_](\d+)[.\-\_]{1,2}(\d+)");
            rst = pattern.findall(line);
            for row in rst:
                chrName, start, stop = row[0], int(row[1]), int(row[2]);
                if fOut != None:
                    print >>fOutH, "%s:%d..%d" % (chrName,start,stop);
                if chrName not in locDict:
                    locDict[chrName] = dict();
                if start in locDict[chrName]:
                    if stop > locDict[chrName][start]:
                        locDict[chrName][start] = stop;
                else:
                    locDict[chrName][start] = stop;
            if re.match("\*", line):
                break;
    elif option1 == 1:
        tmpStr = '';
        for line in fInH:
            if re.match("\*", line):
                break;
            tmpStr += line;
        locInfoDictF, locInfoDictR = dict(), dict();
        if type == "Gene":
            Coln_geneId, Coln_posGene, Coln_posCDS, Coln_posIntron, Coln_posUTR5, \
                Coln_posUTR3 = 0, 1, 2, 3, 4, 5;
            paraLst = tmpStr.split(">");
            for para in paraLst:
                if para == "":
                    continue;
                lineLst = para.splitlines();
                geneFamily = lineLst[0];
                locInfoDictF[geneFamily] = dict();
                for line in lineLst[1:]:
                    if line == "":
                        continue;
                    ary = line.split("\t");
                    assert len(ary) == 6;
                    geneId = ary[Coln_geneId];
                    posGeneLst = ary[Coln_posGene].split(";");
                    posCDSLst = ary[Coln_posCDS].split(";");
                    posIntronLst = ary[Coln_posIntron].split(";");
                    posUTR5Lst = ary[Coln_posUTR5].split(";");
                    posUTR3Lst = ary[Coln_posUTR3].split(";");
                    locInfoDictF[geneFamily][geneId] = {\
                        "Gene":posGeneLst, "CDS":posCDSLst, \
                        "Intron":posIntronLst, "UTR5":posUTR5Lst, \
                        "UTR3":posUTR3Lst };
                    for posGene in posGeneLst:
                        locInfoDictR[posGene] = [geneFamily, geneId, "Gene"];
                    for posCDS in posCDSLst:
                        locInfoDictR[posCDS] = [geneFamily, geneId, "CDS"];
                    for posIntron in posIntronLst:
                        locInfoDictR[posIntron] = [geneFamily, geneId, "Intron"];
                    for posUTR5 in posUTR5Lst:
                        locInfoDictR[posUTR5] = [geneFamily, geneId, "UTR5"];
                    for posUTR3 in posUTR3Lst:
                        locInfoDictR[posUTR3] = [geneFamily, geneId, "UTR3"];
        elif type == "Chromosome" or type == "Probe":
            Coln_pos = 0;
            paraLst = tmpStr.split(">");
            for para in paraLst:
                if para == "":
                    continue;
                lineLst = para.splitlines();
                tag = lineLst[0];
                locInfoDictF[tag] = dict();
                for line in lineLst[1:]:
                    if line == "":
                        continue;
                    ary = line.split("\t");
                    assert len(ary) == 1;
                    posLst = ary[Coln_pos].split(";");
                    for pos in posLst:
                        locInfoDictR[pos] = [tag];
        if option2 == "F":
            locDict = locInfoDictF;
        elif option2 == "R":
            locDict = locInfoDictR;
    return locDict;

def delLocOverlap(locDict):
    lstLoc = locDict.keys();
    locDictNew = locDict;
    flagExit = 1;
    for i in range(0, len(lstLoc)):
        if lstLoc[i] not in locDict:
            continue;
        start1 = locDict[lstLoc[i]][0];
        stop1 = locDict[lstLoc[i]][1];
        if start1 > stop1:
            start1, stop1 = stop1, start1;
        for j in range(i+1, len(lstLoc)):
            if lstLoc[j] not in locDict:
                continue;
            start2 = locDict[lstLoc[j]][0];
            stop2 = locDict[lstLoc[j]][1];
            if start2 > stop2:
                start2, stop2 = stop2, start2;
            if start1 in range(start2,stop2-10) or stop1 in range(start2+10,stop2):
                if lstLoc[i] in locDictNew:
                    del(locDictNew[lstLoc[i]]);
                elif lstLoc[j] in locDictNew:
                    del(locDictNew[lstLoc[j]]);
                flagExit = 0;
                #print start1, stop1, start2, stop2;
    if flagExit == 0:       
        return delLocOverlap(locDictNew);
    else:
        return lstLoc;  

def locCluster(lstLoc):
    flagExit = 1;
    while flagExit == 0:
        for i in range(0, len(lstLoc)):
            locI = sum(lstLoc[i]) / float(len(lstLoc[i]));
            for j in range(i+1, len(lstLoc)):
                locJ = sum(lstLoc[j]) / float(len(lstLoc[j]));
                if abs(locI -locJ) < 2000:
                    locNew = lstLoc[i]+lstLoc[j];
                    lstLoc.pop(i);
                    lstLoc.pop(j-1);
                    lstLoc.append(locNew);
                    #print locI, locJ, locNew;
                    flagExit = 0;
                    continue;
    return lstLoc;

def parsePosStr(posConvStr):
    posLst = [];
    rstAll = re.match("^join\((.+)\)$", posConvStr);
    assert rstAll != None;
    lstPosStr = re.split("; ", rstAll.group(1));
    for posStr in lstPosStr:
        lst = re.match("^(\w+)\:(\d+)..(\d+)$", posStr);
        flagRevComp = 0;
        flagGap = 0;
        if lst == None:
            lst = re.match("^revcom\((\w+)\:(\d+)..(\d+)\)$", posStr);
            flagRevComp = 1;
        assert lst != None;
        if flagGap == 0:
            acc, start, stop = lst.group(1), int(lst.group(2)), int(lst.group(3));
            posLst.append([acc,start,stop]);
    return posLst;

class posConv():
    """#change location on BAC to location on chromosome, e.g.:
    """
    def __init__(self, ):
        pass;
    def BAC2Chr(self, acc, lStart, lStop):
        self.cursor = conn.cursor(MySQLdb.cursors.DictCursor);
        if lStart > lStop:
            lStart, lStop = lStop, lStart;
        queryString = "SELECT * FROM mt3 WHERE acc LIKE '{0}%'";
        self.cursor.execute(queryString.format(acc));
        rst = self.cursor.fetchall();
        if len(rst) == 0:
            return 0;
        elif len(rst) >= 2:
            print "Found multiple entries for '%s'" %acc;
            sys.exit(1);
        row = rst[0];
        bacStart, bacStop, chrStart, chrStop = int(row["bac_start"]), \
            int(row["bac_stop"]), int(row["chr_start"]), int(row["chr_stop"]);
        gStart = self.posConversion(lStart, bacStart, bacStop, \
            chrStart, chrStop, row["orientation"]);
        gStop = self.posConversion(lStop, bacStart, bacStop, \
            chrStart, chrStop, row["orientation"]);
        if gStart > gStop:
            gStart, gStop = gStop, gStart;
        rst = "MtChr%d:%d..%d" % (int(row["chr"]), gStart, gStop);
        return rst;
    def chr2BAC(self, chr, gStart, gStop):
        self.cursor = conn.cursor(MySQLdb.cursors.DictCursor);
        if gStart > gStop:
            gStart, gStop = gStop, gStart;
        queryString = "SELECT * FROM mt3 WHERE chr=%d AND \
            ( (chr_start<=%d AND chr_stop>=%d) \
            OR (chr_start>%d AND chr_stop<%d) \
            OR (chr_start<=%d AND chr_stop>=%d) )" \
            % (chr, gStart, gStart, gStart, gStop, gStop, gStop);
        self.cursor.execute(queryString);
        result = self.cursor.fetchall();
        assert len(result) > 0;
        seqRcdIdLst = [];
        for row in result:
            gStartP, gStopP = gStart, gStop;
            if(row["chr_start"]<=gStart and row["chr_stop"]>=gStart):
                #print "Start:" + row["acc"];
                if(row["chr_start"]<=gStop and row["chr_stop"]>=gStop):
                    pass;
                    #print "Stop:" + row["acc"];
                else:
                    gStopP = row["chr_stop"];
            elif(row["chr_start"]<=gStop and row["chr_stop"]>=gStop):
                gStartP = row["chr_start"];
                #print "Stop" + row["acc"];
            else:
                #print "middle:" + row["acc"];
                gStartP, gStopP = row["chr_start"], row["chr_stop"];
            chrStart, chrStop = int(row["chr_start"]), int(row["chr_stop"]);
            if(row["type"] == "N"):
                gapStart, gapStop = 1, int(row["acc"]);
                lStart = self.posConversion(gStartP, chrStart, chrStop, \
                    gapStart, gapStop, row["orientation"]);
                lStop = self.posConversion(gStopP, chrStart, chrStop, \
                    gapStart, gapStop, row["orientation"]);
                print "\tGap\t(MtChr%d:%d..%d)" % (chr, gStart, gStop);
                sid = "gap(%d)" % (lStop-lStart+1);
            else:
                bacStart, bacStop = int(row["bac_start"]), int(row["bac_stop"]);
                lStart = self.posConversion(gStartP, chrStart, chrStop, \
                    bacStart, bacStop, row["orientation"]);
                lStop = self.posConversion(gStopP, chrStart, chrStop, \
                    bacStart, bacStop, row["orientation"]);
                if lStart > lStop:
                    lStart, lStop = lStop, lStart;
                acc = row["acc"][:row["acc"].index(".")];
                if(row["orientation"] == 1):
                    sid = "%s:%d..%d" % (acc, lStart, lStop);
                elif(row["orientation"] == -1):
                    sid = "revcom(%s:%d..%d)" % (acc, lStart, lStop);
                else:
                    print "Error!";
                    sys.exit(0);
            seqRcdIdLst.append(sid);
        return "join("+"; ".join(seqRcdIdLst)+")";
    def posConversion(self, qry, start1, stop1, start2, stop2, orient):
        if orient==1 or orient==0:
            converted = qry - start1 + start2;
        elif orient == -1:
            converted = stop2 - (qry - start1);
        else:
            print "unknow orientation: %d" % orient;
            sys.exit(0);
        return converted;

if __name__ == "__main__":
    #test for position conversion system
    locDict = my.getLoc(fIn);
    myPosConv = my.posConv();
    fOutH = open(fOut,"w");
    for chrName, value in sorted(locDict.iteritems()):
        for start, stop in sorted(value.iteritems()):
            print "%s:%d..%d\t" % (chrName,start,stop),
            print >>fOutH, "%s:%d..%d\t" % (chrName,start,stop);
            chr = int(chrName[-1]);
            tmp = myPosConv.chr2BAC(chr,start,stop);
            print "%s\t" % tmp;
            for ele in my.parsePosStr(tmp):
                print "\t"+myPosConv.BAC2Chr(ele[0],ele[1],ele[2]),
            print "\n",

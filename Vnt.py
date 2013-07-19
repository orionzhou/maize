# -*- coding: utf-8 -*-
import os.path;

class VntRet():
    def __init__(self, lstAcc=None, delim="\t"):
        self.delim = delim;
        self.Number_acc = 12;
        self.Len_acc = 7;
        self.Coln_id = 1;
        self.Coln_info = 0;
        self.Len_info = 16
        self.Coln_acc = 16;
        self.Coln_sum = self.Coln_acc + self.Number_acc * self.Len_acc;
        self.Len_sum = 8;
        self.filePos1 = 0;
        self.filePos2 = 0;
        self.accColnDict = dict();
        tmp = os.path.join(DIR_Data, "Variant_Mt", "MtChr1_snp");
        tmpH = open(tmp, "r");
        firstLine = tmpH.readline();
        ary = firstLine.strip("\n").split(self.delim);
        if lstAcc == None:
            for cntAcc in range(self.Number_acc):
                acc = ary[self.Coln_acc+cntAcc*self.Len_acc];
                self.accColnDict[acc] = ary.index(acc);
        else:
            for acc in lstAcc:
                assert acc in ary;
                self.accColnDict[acc] = ary.index(acc);
    def retByLine(self, chrName, lineLst, fOut, vntType = "snp"):
        lineLst = sorted(lineLst);
        self.fIn1 = os.path.join(DIR_Data, "Variant_Mt", "%s_%s" % (chrName,vntType));
        self.fInH1 = open(self.fIn1, "r");
        firstLine = self.fInH1.readline();
        cntLine = 0;
        fOutH = open(fOut, "w");
        self.writeHead(fOutH);
        id = "%s:%d..%d" % (chrName,lineLst[0],lineLst[-1]);
        for i in range(lineLst[-1]):
            line = self.fInH1.readline();
            cntLine += 1;
            if cntLine in lineLst:
                lineStripped = line.strip("\n");
                if lineStripped != "":
                    self.writeSNP(lineStripped, fOutH, id);
    def writeHead(self, fOutH):
        print >>fOutH, "\t".join(["VntId", "Class", "Region", "Position", \
            "PosOri", "RefAllele", "VarAllele", "Accession", \
            "NumReadsWithAllele", "Coverage", "Freq", "UniqAlns", "AvgQual"]);
    def writeSNP(self, line, fOutH, id):
        ary = line.split(self.delim);
        if ary[self.Coln_info+2] == "S":
            rst = re.match("^(\d+)$", ary[self.Coln_info+3]); 
        elif ary[self.Coln_info+2] == "I":
            rst = re.match("^(\d+)\^(\d+)$", ary[self.Coln_info+3]);
        elif ary[self.Coln_info+2] == "D":
            rst = re.match("^\[(\d+)\.\.(\d+)\]$", ary[self.Coln_info+3]);
        assert rst != None;
        pos = int(rst.group(1));
        dict_row = { "Reference": ary[self.Coln_info+0],
                    "Variant": ary[self.Coln_info+1],
                    "Class": ary[self.Coln_info+2],
                    "Position": ary[self.Coln_info+3],
                    "RefAllele": ary[self.Coln_info+4],
                    "VarAllele": ary[self.Coln_info+5],
                    "MatchedVariant": ary[self.Coln_sum+4],
                    "MatchedReference": ary[self.Coln_sum+5],
                    "BelowThreshold": ary[self.Coln_sum+6],
                    "NoCoverage": ary[self.Coln_sum+7] 
            };
        if id == "":
            id = dict_row["Reference"];
        for acc in sorted(self.accColnDict.keys()):
            accColn = self.accColnDict[acc];
            accInfo = { "Summary": ary[accColn], \
                        "NumReadsWithAllele": ary[accColn+1], \
                        "Coverage": ary[accColn+2], \
                        "Freq": ary[accColn+3], \
                        "UniqAlns": ary[accColn+4], \
                        "AvgQual": ary[accColn+5], \
                        "MaxQual": ary[accColn+6], \
                };
            print >>fOutH, "\t".join([dict_row["Variant"], dict_row["Class"], \
                id, str(pos), dict_row["Position"], dict_row["RefAllele"], \
                dict_row["VarAllele"], acc[0:5], accInfo["NumReadsWithAllele"], accInfo["Coverage"], \
                accInfo["Freq"], accInfo["UniqAlns"], accInfo["AvgQual"]]);
    def retByLoc(self, locDict, fOut, vntType="all"):
        #resultDict = dict();
        fOutH = open(fOut, "w");
        self.writeHead(fOutH);
        for chrName, value in sorted(locDict.iteritems()):
            print "%s:" % chrName;
            lstLocByChr = [];
            for start, stop in sorted(value.iteritems()):
                lstLocByChr.append([chrName,start,stop]);
            #resultDict.update(self.retByChr(lstLocByChr, chrName));
            self.retByChr(lstLocByChr, chrName, fOutH, vntType);
        #with open(fOut, "w") as fOutH:
            #pickle.dump(resultDict, fOutH);
    def retByChr(self, lstLocByChr, chrName, fOutH, vntType):
        self.fIn1 = os.path.join(DIR_Data, "Variant_Mt", "%s_snp" % chrName);
        self.fIn2 = os.path.join(DIR_Data, "Variant_Mt", "%s_indel" % chrName);
        self.fInH1 = open(self.fIn1, "r");
        self.fInH2 = open(self.fIn2, "r");
        firstLine1 = self.fInH1.readline();
        firstLine2 = self.fInH2.readline();
        assert firstLine1 == firstLine2;
        self.filePos1 = self.fInH1.tell();
        self.filePos2 = self.fInH2.tell();
        cnt = 0;
        start_prev = 0;
        dictRst = dict();
        for ele in lstLocByChr:
            cnt += 1;
            chrName2, start, stop = ele[0], ele[1], ele[2];
            assert chrName == chrName2 and start_prev <= start;
            id = "%s:%d..%d" % (chrName,start,stop);
            print "\t%s" % id,
            if vntType == "snp" or vntType == "all":
                self.retOne(chrName, start, stop, 1, id, fOutH);
            if vntType == "indel" or vntType == "all":
                self.retOne(chrName, start, stop, 2, id, fOutH);
            print "\n",
            start_prev = start; 
    def retOne(self, chrName, start, stop, option, id, fOutH):
        if option == 1:
            fInH = self.fInH1;
            print "\tsnp: ",
        elif option == 2:
            fInH = self.fInH2;
            print "\tindel: ",
        pos_prev, pos_curr, cntVnt, cntPos= 0, 0, 0 ,0;
        dict_parsed = dict();
        for line in fInH:
            line = line.strip("\n");
            if line == "":
                break;
            ary = line.split(self.delim);
            if ary[0] != chrName:
                print line;
            assert ary[0] == chrName;
            if ary[self.Coln_info+2] == "S":
                rst = re.match("^(\d+)$", ary[self.Coln_info+3]); 
            elif ary[self.Coln_info+2] == "I":
                rst = re.match("^(\d+)\^(\d+)$", ary[self.Coln_info+3]);
            elif ary[self.Coln_info+2] == "D":
                rst = re.match("^\[(\d+)\.\.(\d+)\]$", ary[self.Coln_info+3]);
            assert rst != None;
            pos_curr = int(rst.group(1));
            if pos_curr < start:
                continue;
            elif pos_curr>=start and pos_curr<=stop:
                cntVnt += 1;
                if pos_curr != pos_prev:
                    cntPos += 1;
                self.writeSNP(line, fOutH, id);
                pos_prev = pos_curr;
            elif pos_curr > stop:
                if option == 1:
                    #print "(%d-" % self.fInH1.tell(),
                    curPos = self.fInH1.tell();
                    if curPos-len(line)-90000 < 0:
                        self.fInH1.seek(0, 0);
                    else:
                        self.fInH1.seek(-len(line)-90000, 1);
                    tmp = self.fInH1.readline();
                    self.filePos1 = self.fInH1.tell();
                    #print "%d)" % self.filePos1,
                elif option == 2:
                    curPos = self.fInH2.tell();
                    if curPos-len(line)-90000 < 0:
                        self.fInH2.seek(0, 0);
                    else:
                    #print "(%d-" % self.fInH2.tell(),
                        self.fInH2.seek(-len(line)-90000, 1);
                    tmp = self.fInH2.readline();
                    self.filePos2 = self.fInH2.tell();
                    #print "%d)" % self.filePos2,
                break;
        print "{0}({1} UniqPos)".format(cntVnt, cntPos),

class VntReport():
    def __init__(self, cutOffCov=2, cutOffFreq=0.7, cutOffUniq=2, mode="chr", fLocInfo=None):
        self.cutOffCov = cutOffCov;
        self.cutOffFreq = cutOffFreq;
        self.cutOffUniq = cutOffUniq;
        assert mode in ["gene", "chr"];
        self.mode = mode;
        self.locInfoDict = None;
        self.vntEffectDict = dict();
        if self.mode == "gene":
            assert fLocInfo != None;
            self.locInfoDict = getLoc(fLocInfo, 1, "R", 1);
    def getAcc(self, fVnt):
        lstAcc = [];
        with open(fVnt, "r") as fVntH:
            firstLine = fVntH.readline();
            firstPos = 0;
            firstVnt = "";
            cnt = 0;
            for line in fVntH:
                ary = line.strip("\n").split("\t");
                if cnt == 0:
                    firstPos = int(ary[3]);
                    firstVnt = ary[6];
                    lstAcc.append(ary[7]);
                else:
                    if int(ary[3]) == firstPos and ary[6] == firstVnt:
                        lstAcc.append(ary[7]);
                    else:
                        break;
                cnt += 1;
            assert cnt == len(lstAcc);
        print "%d accessions:" % len(lstAcc);
        print "\t"+" ".join(lstAcc);
        self.lstAcc = lstAcc;
    def readIn(self, fVnt, fReport, fReport2="None"):
#       if self.mode == "gene":
#           assert fReport2 != None;
#           self.fRH2 = open(fReport2, "w");
        self.fRH = open(fReport, "w");
        self.writeHead();
        self.getAcc(fVnt);
        accDict, cntCovDict, cntVntDict, cntVntCalledDict, cntVntCoveredDict = \
            dict(), dict(), dict(), dict(), dict();
        for acc in self.lstAcc:
            accDict[acc], cntCovDict[acc], cntVntDict[acc], \
                cntVntCalledDict[acc], cntVntCoveredDict[acc] = "",0,0,0,0;
        fVntH = open(fVnt, "r");
        firstLine = fVntH.readline();
        locLst, cntLine, cntVnt, cntLoc =[], 0, 0, 0;
        cntAcc = 0;
        regionPrev = "";
        for line in fVntH:
            line = line.strip("\n");
            if line == "":
                break;
            ary = line.split("\t");
            vntRow = {  "VntId":ary[0], "Class":ary[1], \
                        "Region":ary[2],
                        "Position": int(ary[3]), \
                        "PosOri": ary[4],
                        "RefAllele": ary[5], \
                        "VarAllele": ary[6], \
                        "Accession": ary[7], \
                        "NumReadsWithAllele": int(ary[8]), \
                        "Coverage": int(ary[9]), \
                        "Freq": float(ary[10]), \
                        "UniqAlns": int(ary[11])
                };
            if ary[12] == "":
                vntRow["AvgQual"] = 0;
            else:
                vntRow["AvgQual"] = float(ary[12]);
            acc = vntRow["Accession"];
            assert acc in self.lstAcc;
            if regionPrev != vntRow["Region"] and regionPrev != "":
                cntVnt = cntLine / len(self.lstAcc);
                self.locLst, self.accDict, self.cntVnt, self.cntLoc = \
                    locLst, accDict, cntVnt, cntLoc;
                self.writeStat(regionPrev, cntCovDict, cntVntDict, \
                    cntVntCalledDict, cntVntCoveredDict);
                locLst, cntLine, cntVnt, cntLoc =[], 0, 0, 0;
                for acc2 in self.lstAcc:
                    accDict[acc2], cntCovDict[acc2], cntVntDict[acc2], \
                        cntVntCalledDict[acc2], cntVntCoveredDict[acc2]= "", 0, 0, 0, 0;
            regionPrev = vntRow["Region"];
            cntLine += 1;
            if vntRow["Position"] not in locLst:
                #print vntRow["Position"], len(locLst);
                #for a,b in accDict.iteritems():
                    #assert len(locLst) == len(b);
                locLst.append(vntRow["Position"]);
                cntLoc += 1;
                cntAcc = 0;
            cntAcc += 1;
            #be careful of multi-variants at the same position
            if cntAcc <= len(self.lstAcc): 
                cntCovDict[acc] += vntRow["Coverage"];
                cntVntDict[acc] += 1;
                if vntRow["Coverage"] >= self.cutOffCov:
                    #assert len(vntRow["RefAllele"])==1 and len(vntRow["VarAllele"])==1;
                    cntVntCoveredDict[acc] += 1;
                    if vntRow["UniqAlns"]>=self.cutOffUniq and \
                            vntRow["Freq"]>=self.cutOffFreq:
                        accDict[acc] += vntRow["VarAllele"].upper();
                        cntVntCalledDict[acc] += 1;
#                       if self.mode == "gene": self.vntEffect(vntRow);
                    else:
                        accDict[acc] += vntRow["RefAllele"].upper();
                else:
                    accDict[acc] += "N";
            else:
                if vntRow["Coverage"]>=self.cutOffCov and \
                        vntRow["UniqAlns"]>=self.cutOffUniq \
                        and vntRow["Freq"]>=self.cutOffFreq:
                    accDict[acc] = accDict[acc][0:-1] + vntRow["VarAllele"].upper();
                    #accDict[acc] = accDict[acc][0:-1] + vntRow["RefAllele"].upper(); 
                    cntVntCalledDict[acc] += 1;
#                   if self.mode == "gene": self.vntEffect(vntRow);
        #deal with the last region
        cntVnt = cntLine / len(self.lstAcc);
        self.locLst, self.accDict, self.cntVnt, self.cntLoc = \
            locLst, accDict, cntVnt, cntLoc;
        self.writeStat(regionPrev, cntCovDict, cntVntDict, cntVntCalledDict, \
            cntVntCoveredDict);
        fVntH.close();
        self.fRH.close();
    def writeHead(self):
        if self.mode == "gene":
            print >>self.fRH, "\t".join(["geneFamily", "geneId", \
                "type", "region", "acc", \
                "totCov", "totVnt", "calledVnt", "coveredVnt", "length"]);
            print >>self.fRH2, "\t".join(["geneFamily", "geneId", \
                "type", "region", "acc", "position", "class", "desc", \
                "effect"]);
        else:
            print >>self.fRH, "\t".join(["chr", "posMiddle", \
                "region", "acc", \
                "totCov", "totVnt", "calledVnt", "coveredVnt", "length"]);
    def writeStat(self, region, cntCovDict, cntVntDict, cntVntCalledDict, \
        cntVntCoveredDict):
        print "{0}: {1} variants / {2} unique loci".format(region, \
            self.cntVnt, self.cntLoc);
        assert self.cntLoc == len(self.locLst);
        #for acc,hap in self.accDict.iteritems():
            #print acc, len(hap), self.cntLoc;
            #assert self.cntLoc == len(hap);
        tmp = region.split(":");
        tmp2 = tmp[1].split("..");
        assert len(tmp2) == 2;
        length = int(tmp2[1])-int(tmp2[0])+1;
        if self.mode == "gene":
            assert region in self.locInfoDict;
            locInfo = self.locInfoDict[region];
            geneFamily, geneId, type = locInfo[0], locInfo[1], locInfo[2];
            for acc in self.lstAcc:
                print >>self.fRH, "\t".join([geneFamily, geneId, type, region, acc, \
                    str(cntCovDict[acc]), str(cntVntDict[acc]), \
                    str(cntVntCalledDict[acc]), str(cntVntCoveredDict[acc]), \
                    str(length)]);
        else:
            chr, posMiddle = tmp[0], (int(tmp2[0])+int(tmp2[1])) / 2;
            for acc in self.lstAcc:
                print >>self.fRH, "\t".join([chr, str(posMiddle), region, acc, \
                    str(cntCovDict[acc]), str(cntVntDict[acc]), \
                    str(cntVntCalledDict[acc]), str(cntVntCoveredDict[acc]), \
                    str(length)]);
    def vntEffect(self, vntRow):
        region = vntRow["Region"]
        assert region in self.locInfoDict;
        locInfo = self.locInfoDict[region];
        geneFamily, geneId, type = locInfo[0], locInfo[1], locInfo[2];
        tmp = region.split(":");
        chrName = tmp[0];
        vntKey = "%s:%d:%s:%s" % (chrName, vntRow["Position"], \
            vntRow["RefAllele"], vntRow["VarAllele"]);
        if vntKey in self.vntEffectDict:
            [effect, desc] = self.vntEffectDict[vntKey];
        else:
            myDesc = VntDesc();
            [effect, desc] = myDesc.desc(vntKey, type);
            self.vntEffectDict[vntKey] = [effect, desc];
        print >>self.fRH2, "\t".join([geneFamily, geneId, type, region, \
            vntRow["Accession"], str(vntRow["Position"]), vntRow["Class"], \
            desc, effect]);
    def write(self, cutOff_MinorAlleleCount, fOutLst):
        print "---->output";
        tag = [os.path.exists(fOut) for fOut in fOutLst];
        fOutHLst = [open(fOut,"w") for fOut in fOutLst];
        numLoci = len(self.locLst);
        dictEncoding1 = {"A":"1", "C":"2", "G":"3", "T":"4", "N":"0"};
        dictEncoding2 = {"A":"A", "C":"C", "G":"G", "T":"T", "N":"?"};
        dictEncoding3 = {"A":"L", "C":"O", "G":"P", "T":"Q", "N":"Z"};
        pAlleleDict, pLocLst = dict(), [];
        for acc in self.lstAcc:
            pAlleleDict[acc] = "";
        cntMono, cntBi, cntMulti, cntP = 0, 0, 0, 0;
        for i in range(numLoci):
            tmpDict = dict();
            for acc, hap in sorted(self.accDict.iteritems()):
                tmpDict[acc] = hap[i];
            alleleCntDict = dict();
            for allele in tmpDict.values():
                if allele != "N":
                    if allele not in alleleCntDict:
                        alleleCntDict[allele] = 0;
                    alleleCntDict[allele] += 1;
            n = len(alleleCntDict.keys());
            if n == 1:
                cntMono += 1;
            elif n == 2:
                cntBi += 1;
                flag = True;
                for allele,alleleCnt in alleleCntDict.iteritems():
                    if alleleCnt < cutOff_MinorAlleleCount:
                        flag = False;
                if flag == True:
                    cntP += 1;
                    print >>fOutHLst[1], "L%04d\t%d" % (cntBi, self.locLst[i]);
                    pLocLst.append(self.locLst[i]);
                    for acc in self.lstAcc:
                        pAlleleDict[acc] += tmpDict[acc];
            elif n > 2:
                print " Position(%d) has more than 2 alleles -> discarded" \
                    % self.locLst[i];
                cntMulti += 1;
        print >>fOutHLst[2], len(self.lstAcc);
        print >>fOutHLst[2], cntP;
        print >>fOutHLst[2], "S" * len(pLocLst);
        print >>fOutHLst[2], " ".join([str(loc) for loc in pLocLst]);
        locMin, locMax = min(pLocLst), max(pLocLst);
        assert locMin == pLocLst[0] and locMax == pLocLst[-1];
        print >>fOutHLst[5], "%d\t%d\tL" % \
            (cntP, round((locMax - locMin + 1)/float(1000)));
        print >>fOutHLst[5], \
            "\n".join(["%.3f" % (int(loc-locMin+1)/float(1000)) for loc in pLocLst]);
        cntIndi = 0;
        print >>fOutHLst[4], " %d %d 1" % (len(self.lstAcc), cntP);
        for acc in self.lstAcc:
            assert cntP == len(pAlleleDict[acc]);
            cntIndi += 1;
            tmpStr1 = "%s\t%s\t0\t0\t2\t0\t{0}" % (cntIndi, acc);
            strConverted1, strConverted2 = "", "";
            for ch in pAlleleDict[acc]:
                assert ch in dictEncoding1;
                strConverted1 += "{0} {0}  ".format(dictEncoding1[ch]);
                strConverted2 += dictEncoding2[ch];
            print >>fOutHLst[0], tmpStr1.format(strConverted1);
            print >>fOutHLst[2], acc;
            print >>fOutHLst[2], strConverted2;
            print >>fOutHLst[2], strConverted2;
            print >>fOutHLst[3], ">" + acc;
            print >>fOutHLst[3], pAlleleDict[acc], "\n";
            print >>fOutHLst[4], ">" + acc;
            print >>fOutHLst[4], strConverted2;
        for i in range(len(tag)):
            fOutHLst[i].close();
            if tag[i] == False:
                os.chmod(fOutLst[i], 0777);
        print "---->summary";
        print " {0} Variants / {1} UniqLoci: ".format(self.cntVnt, self.cntLoc);
        print " {0}(mono-allelic) + {1}(bi-allelic) + {2}(multi-allelic)".format(\
            cntMono, cntBi, cntMulti);
        print " {0} picked with cutOff(Minor Allele Count) = {1}".\
            format(cntP, cutOff_MinorAlleleCount);
    def runRsq(self, fIn, fOut, fTmp="rsq.tmp"):
        os.system("rsq -i {0} -c 4 >> {1}".format(fIn, fTmp));
        fOutH = open(fOut, "w");
        with open(fTmp, "r") as fTmpH:
            for line in fTmpH:
                line = line.strip("\n");
                if line[0] == "#":
                    print >>fOutH, line.strip("\n");
                elif line[0] == "p":
                    print >>fOutH, "\t".join(["posi","posj","distance","rsq","D","Dprime"]);
                elif line != "":
                    ary = line.split("\t");
                    assert len(ary) == 5;
                    pos1,pos2 = int(ary[0]), int(ary[1]);
                    distance = abs(self.locLst[pos1-1]-self.locLst[pos2-1]);
                    print >>fOutH, "\t".join([ary[0],ary[1],str(distance),ary[2],ary[3],ary[4]]);
        os.remove(fTmp);

class VntDesc():
    def __init__(self):
        self.cursor = conn.cursor(MySQLdb.cursors.DictCursor);
    def desc(self, vntKey, type):
        tmp = vntKey.split(":");
        assert len(tmp) == 4;
        chrName, posG, ref, var = tmp[0], int(tmp[1]), tmp[2], tmp[3];
        if type == "Gene":
            qryString = "SELECT ID,Start,Stop FROM mt_hapmap.mt3_annotation \
                WHERE SeqID='{0}' AND Start<={1} AND Stop>={1} AND Type='{2}'";
            self.cursor.execute(qryString.format(chrName, posG, "CDS"));
            rst = self.cursor.fetchall();
            if len(rst) >= 1:
                type = "CDS";
            else:
                qryString = "SELECT ID,Start,Stop FROM mt_hapmap.mt3_annotation \
                    WHERE SeqID='{0}' AND Start<={1} AND Stop>={1} AND Type='{2}'";
                self.cursor.execute(qryString.format(chrName, posG, "exon"));
                rst = self.cursor.fetchall();
                if len(rst) >= 1:
                    type = "UTR5";
                else:
                    type = "Intron";
        assert type in ["CDS", "Intron", "UTR5", "UTR3"];
        if type == "Intron" or type == "UTR5" or type == "UTR3":
            effect = "-";
            desc = "%s/%s" % (ref, var);
        elif type == "CDS":
            [effect, desc] = self.descCDS(chrName, posG, ref, var);
            self.cursor.execute(qryString.format(chrName, posG, "Exon"));
        return [effect, desc];
    def descCDS(self, chrName, posG, ref, var):
        qryString = "SELECT Parent,Start,Stop,Strand FROM mt_hapmap.mt3_annotation \
            WHERE SeqID='{0}' AND Start<={1} AND Stop>={1} AND Type='{2}'";
        self.cursor.execute(qryString.format(chrName, posG, "CDS"));
        rst = self.cursor.fetchall();
        assert len(rst) >= 1;
        row = rst[0];
        orient, mRNAId = row["Strand"], row["Parent"];
        qryString = "SELECT ID,SeqID,Start,Stop,Strand FROM \
                mt_hapmap.mt3_annotation WHERE \
                Parent='{0}' AND Type='{1}'";
        self.cursor.execute(qryString.format(mRNAId, "CDS"));
        rst = self.cursor.fetchall();
        assert len(rst) >= 1;
        posDict = dict();
        myret = SeqRet();
        cdsStr = "";
        posR, totLen= 0, 0;
        flagHit = False;
        for row in rst:
            assert row["Strand"] == orient;
            posDict[row["Start"]] = row["Stop"];
        if orient == "+":
            for start, stop in sorted(posDict.iteritems()):
                seqRcd = myret.ret(int(chrName[-1]), start, stop);
                cdsStr += seqRcd.seq.tostring();
                if posG>=start and posG<=stop:
                    posR = posG-start+1 + totLen;
                    flagHit = True;
                totLen += (stop-start+1);
        elif orient == "-":
            for start, stop in sorted(posDict.iteritems(), reverse=True):
                seqRcd = myret.ret(int(chrName[-1]), start, stop);
                cdsStr += seqRcd.seq.reverse_complement().tostring();
                if posG>=start and posG<=stop:
                    posR = stop-posG+1 + totLen;
                    flagHit = True;
                totLen += (stop-start+1);
        if orient == "-":
            ref = Seq(ref, IUPAC.unambiguous_dna).reverse_complement().tostring();
            var = Seq(var, IUPAC.unambiguous_dna).reverse_complement().tostring();
        if totLen % 3 != 0 or flagHit != True:
            print vntKey, posR, "/", totLen, ref, cdsStr[posR-1];
        assert totLen % 3 == 0 and flagHit == True;
        cntAA = totLen / 3;
        posR_AA, phase = int((posR-1)/3), (posR-1)%3;
        codonRef = Seq(cdsStr[posR_AA*3:posR_AA*3+3], IUPAC.unambiguous_dna);
        codonVar = Seq(cdsStr[posR_AA*3:posR-1]+var+cdsStr[posR:posR_AA*3+3], \
            IUPAC.unambiguous_dna);
        if "N" in codonRef:
            aaRef, aaVar = "x", "x";
        else:
            try:
                aaRef, aaVar= codonRef.translate(table=1).tostring(), \
                    codonVar.translate(table=1).tostring();
            except:
                print vntKey, posR, "/", totLen, ref, cdsStr[posR-1], codonRef;
                print cdsStr;
        desc = "%s->%s|%s->%s" % (codonRef, codonVar, aaRef, aaVar);
        effect = None;
        if aaRef == aaVar:
            effect = "S";
        else:
            effect = "N";
            if aaRef == "*" and aaVar != "*":
                print "\tStop Codon (%s->%s)" % (aaRef, aaVar);
                assert posR_AA == cntAA-1;
                effect += "E";
            elif aaRef != "*" and aaVar == "*":
                print "\tPremature Stop Codon (%s->%s)" % (aaRef, aaVar);
                effect += "P";
            elif posR_AA == 0:
                print "\tStart Codon (%s->%s)" % (aaRef, aaVar);
                assert aaRef == "M" and aaVar != "M";
                effect += "B";
        return [effect, desc];

lstAcc = ["HM001","HM002","HM003","HM004","HM005_gsnapout8mm", \
    "HM006","HM009","HM011","HM015","HM017-I", \
    'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];

cutOffCov, cutOffFreq, cutOffUniq= 2, 0.7, 2;
cutOffMAC = 1;

suffix = "MtChr";
fIn = os.path.join(my.DIR_In, "location").format(suffix);
assert os.path.exists(fIn);
DIR_Work = os.path.join(DIR_Misc, suffix);
if not os.path.exists(DIR_Work):
    os.mkdir(DIR_Work);
fVnt = os.path.join(DIR_Work, "vnt.txt");
fReport = os.path.join(DIR_Work, "vntStat.txt");
fOut1 = os.path.join(DIR_Work, "in_haploview.txt");
fOut2 = os.path.join(DIR_Work, "in_haploview_loc.txt");
fOut3 = os.path.join(DIR_Work, "in_phase.txt");
fOut4 = os.path.join(DIR_Work, "in_rsq.txt");
#fOut42 = os.path.join(DIR_Work, "rsq.out");
fOut5 = os.path.join(DIR_Work, "in_LDhat.txt");
fOut6 = os.path.join(DIR_Work, "in_LDhat_loc.txt");
fOut7 = os.path.join(DIR_Work, "in_LDhot.txt");
fOutLst = [fOut1, fOut2, fOut3, fOut4, fOut5, fOut6, fOut7];

if __name__ == "__main__":   
    locDict = my.getLoc(fIn, 2);
    myVntRet = my.VntRet(lstAcc);
    myVntRet.retByLoc(locDict, fVnt, "snp");
#   myVntOut = my.VntReport(cutOffCov, cutOffFreq, cutOffUniq);
#   myVntOut.readIn(fVnt, fReport);
#   myVntOut.write(cutOffMAC, fOutLst);
#   myVntOut.runRsq(fOut4,fOut42);

#    myDesc = VntDesc();
#    testAry = ["MtChr1:2144758:A:g", "MtChr1:2144754:T:a", "MtChr5:30092376:C:t", "MtChr5:30192504:C:t", "MtChr6:8747283:G:t"];
#    for test in testAry:
#        print test+"\t"+"\t".join(myDesc.desc(test, "Gene"));

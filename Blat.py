import os
import os.path

class blat():
    def __init__(self, qrySeqRcdLst, dbOpt=1, oocOpt=1, idenT=0.9):
        self.qrySeqRcdLst = qrySeqRcdLst;
        self.seqLst, self.lstLen, self.seqDescLst= [], [], [];
        for seqRcd in self.qrySeqRcdLst:
            self.seqLst.append(seqRcd.id);
            self.seqDescLst.append(seqRcd.description);
            self.lstLen.append(len(seqRcd));
        self.blatCmd = "{0}/blat {1} {2} -t={3} -q={4} {5} -noHead {6}";
        assert dbOpt in range(1,5) and oocOpt in range(1,4);
        dbname, oocCmd = "", "";
        if dbOpt == 1:
            dbname = "Mt3.0.fa";
        elif dbOpt == 2:
            dbname = "Mt3.0_chr1-8.fa";
        elif dbOpt == 3:
            dbname = "Mt3.0_cds.fa";
        elif dbOpt == 4:
            dbname = "Mt3.0_proteins.fa";
        if oocOpt == 1:
            oocCmd = "-ooc=%s" % os.path.join(DIR_Data, "Genome_Mt", "11.ooc");
        elif oocOpt == 3:
            oocCmd = "-makeOoc=11.ooc";
        self.fDb = os.path.join(DIR_Data, "Genome_Mt", dbname);
        self.dbType, self.oocCmd, self.identyThresh = "dna", oocCmd, idenT;
    def blat(self, fQry, fQryType, fOut):
        with open(fQry, "w") as fSeqH:
            SeqIO.write(self.qrySeqRcdLst, fSeqH, "fasta");
        os.system(self.blatCmd.format(DIR_Blat, self.fDb, fQry, self.dbType, \
            fQryType, self.oocCmd, fOut));
    def blatFilter(self, fIn, fOut):
        fInH, fOutH = open(fIn, "r"), open(fOut, "w");
        for line in fInH:
            line = line.strip("\n");
            if line != "":
                ary = line.split("\t");
                if len(ary) != 21: print len(ary), ary;
                assert len(ary) == 21;
                matches, misMatches, repMatches, nCount, qNumbInsert, qBaseInsert, \
                tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, \
                tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts = \
                    int(ary[0]), int(ary[1]), int(ary[2]), int(ary[3]), int(ary[4]), \
                    int(ary[5]), int(ary[6]), int(ary[7]), ary[8], ary[9], \
                    int(ary[10]), int(ary[11]), int(ary[12]), ary[13], \
                    int(ary[14]), int(ary[15]), int(ary[16]), int(ary[17]), \
                    ary[18], ary[19], ary[20];
                blockSizeAry = blockSizes.strip(",").split(",");
                qStartAry = qStarts.strip(",").split(",");
                tStartAry = tStarts.strip(",").split(",");
                assert strand in ["+", "-"];
                assert len(blockSizeAry)==blockCount and len(qStartAry)==blockCount \
                    and len(tStartAry)==blockCount;
                qLocAry, tLocAry, cntMA = [], [], matches;
                x1, x3 = int(blockSizeAry[0]), int(tStartAry[0]);
                tStartG, tEndG = x3+1, x3+x1;
                for i in range(blockCount):
                    x1, x2, x3 = int(blockSizeAry[i]), int(qStartAry[i]), \
                        int(tStartAry[i]);
                    tStartTmp, tEndTmp = x3+1, x3+x1;
                    if tStartTmp < tStartG:
                        tStartG = tStartTmp;
                    if tEndTmp > tEndG:
                        tEndG = tEndTmp;
                    if x1 < 12:
                        cntMA -= x1;
                percentId = int((cntMA / float(qSize)) * 100);
                tRange = "%s:%d..%d" % (tName, tStartG, tEndG);
                tLength = tEndG - tStartG + 1;
            if cntMA/float(qSize) >= self.identyThresh and tLength/float(cntMA) <= 5:
                print >>fOutH, line;
        fInH.close() and fOutH.close();
    def blatPretty(self, fIn, fSeq, fOut):
        cmd = "{0}/pslPretty {1} {2} {3} {4}";
        os.system(cmd.format(DIR_Blat, fIn, self.fDb, fSeq, fOut));
    def blatOut(self, fIn, fOut, fSeqHit, fSeqNoHit):
        fInH, fOutH = open(fIn, "r"), open(fOut, "w");
        cntHitDict = dict();
        for line in fInH:
            line = line.strip("\n");
            if line != "":
                ary = line.split("\t");
                assert len(ary) == 21;
                matches, misMatches, repMatches, nCount, qNumbInsert, qBaseInsert, \
                tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, \
                tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts = \
                    int(ary[0]), int(ary[1]), int(ary[2]), int(ary[3]), int(ary[4]), \
                    int(ary[5]), int(ary[6]), int(ary[7]), ary[8], ary[9], \
                    int(ary[10]), int(ary[11]), int(ary[12]), ary[13], \
                    int(ary[14]), int(ary[15]), int(ary[16]), int(ary[17]), \
                    ary[18], ary[19], ary[20];
                qryIndex = self.seqLst.index(qName);
                qryLen, qryDesc = self.lstLen[qryIndex], self.seqDescLst[qryIndex];
                assert qryLen == qSize;
                if qName not in cntHitDict:
                    cntHitDict[qName] = 0;
                    print >>fOutH, "%s\t%s" % (qName, qryDesc);
                cntHitDict[qName] += 1;
                blockSizeAry = blockSizes.strip(",").split(",");
                qStartAry = qStarts.strip(",").split(",");
                tStartAry = tStarts.strip(",").split(",");
                assert strand in ["+", "-"];
                assert len(blockSizeAry)==blockCount and len(qStartAry)==blockCount \
                    and len(tStartAry)==blockCount;
                qLocAry, tLocAry, cntMA = [], [], matches;
                x1, x3 = int(blockSizeAry[0]), int(tStartAry[0]);
                tStartG, tEndG = x3+1, x3+x1;
                for i in range(blockCount):
                    x1, x2, x3 = int(blockSizeAry[i]), int(qStartAry[i]), \
                        int(tStartAry[i]);
                    tStartTmp, tEndTmp = x3+1, x3+x1;
                    if tStartTmp < tStartG:
                        tStartG = tStartTmp;
                    if tEndTmp > tEndG:
                        tEndG = tEndTmp;
                    if x1 >= 12:
                        tLocAry.append("%s:%d..%d" % (tName, x3+1, x3+x1));
                        if strand == "+":
                            qLocAry.append("Query:%d..%d" % (x2+1, x2+x1));
                        else:
                            qLocAry.append("Query:%d..%d" % (qSize-x2, qSize-x2-x1+1));
                    else:
                        cntMA -= x1;
                percentId = int((cntMA / float(qSize)) * 100);
                tRange = "%s:%d..%d" % (tName, tStartG, tEndG);
                tLength = tEndG - tStartG + 1;
                print >>fOutH, "\t{0}%\t{1}/{2}\t{3}\t{4}\t{5}\t{6}"\
                    .format(percentId, cntMA, qSize, ";".join(tLocAry), \
                    ";".join(qLocAry), tRange, tLength);
        fInH.close() and fOutH.close();
        self.blatSum(cntHitDict, fSeqHit, fSeqNoHit);
    def blatSum(self, cntHitDict, fSeqHit, fSeqNoHit):
        seqRcdHitLst, seqRcdNoHitLst = [], [];
        for seqRcd in self.qrySeqRcdLst:
            if seqRcd.id in cntHitDict:
                seqRcdHitLst.append(seqRcd);
            else:
                seqRcdNoHitLst.append(seqRcd);
        with open(fSeqHit, "w") as fSeqHitH:
            SeqIO.write(seqRcdHitLst, fSeqHitH, "fasta");
        with open(fSeqNoHit, "w") as fSeqNoHitH:
            SeqIO.write(seqRcdNoHitLst, fSeqNoHitH, "fasta");
        assert len(cntHitDict.keys()) == len(seqRcdHitLst);
        print "\t".join(["Hits on Mt3.0", "Count"]);
        print "\t".join(["0", str(len(seqRcdNoHitLst))]);
        hitCntSet = set(cntHitDict.values());
        for hitCnt in hitCntSet:
            print "%s\t%s" % (hitCnt, cntHitDict.values().count(hitCnt));       


suffix = "cds306";
fQry = os.path.join(my.DIR_In, "in_{0}.fa".format(suffix));
assert os.path.exists(fQry);
DIR_Work = os.path.join(my.DIR_Misc, suffix);
if not os.path.exists(DIR_Work):
    os.mkdir(DIR_Work);
fSeq = os.path.join(DIR_Work, "seq.fa");
fBlatRst = os.path.join(DIR_Work, "blat_1.psl");
fBlatFiltered = os.path.join(DIR_Work, "blat_2.psl");
fBlatOut = os.path.join(DIR_Work, "blat_3.txt");
fBlatPretty = os.path.join(DIR_Work, "blat_4.psl");
fSeqHit = os.path.join(DIR_Work, "seqHit.fa");
fSeqNoHit = os.path.join(DIR_Work, "seqNoHit.fa");

if __name__ == '__main__':
    option = 1; #set option=2 for probeSet sequence
    seqRcdLst = my.getBlastQry(fQry, option); 
    dbOption, oocOption, qryType = 1, 1, "dna";
    identyThresh = 0.7;
    myBlat = my.blat(seqRcdLst, dbOption, oocOption, identyThresh);
    #myBlat.blat(fSeq, qryType, fBlatRst);
    #myBlat.blatFilter(fBlatRst, fBlatFiltered);
    #myBlat.blatOut(fBlatFiltered, fBlatOut, fSeqHit, fSeqNoHit);
    #myBlat.blatPretty(fBlatFiltered, fSeq, fBlatPretty);

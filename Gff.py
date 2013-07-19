import re
import os

DIR_Genome
fInA1 = os.path.join(DIR_Annotation, "Mt3.0_a_chr.gff");
fInA2 = os.path.join(DIR_Annotation, "Mt3.0_a_pgp.txt");
fOutA = os.path.join(DIR_Annotation, "Mt3.0_a_BAC_o.gff");

fInB = os.path.join(DIR_Annotation, "Mt3.0_a_BAC_info.txt");
fOutB = os.path.join(DIR_Annotation, "Mt3.0_a_BAC_u.gff");

fOutC = os.path.join(DIR_Annotation, "Mt3.0_assembly.gff");

fInD = os.path.join(DIR_Annotation, "medicago3.0_20090702_MIPSGFF.gff");
fOutD1 = os.path.join(DIR_Annotation, "Mt3.01.gff");
fOutD2 = os.path.join(DIR_Annotation, "Mt3.02.gff");
fOutD3 = os.path.join(DIR_Annotation, "Mt3.03_phase.txt");
fOutD4 = os.path.join(DIR_Annotation, "Mt3.04.gff");
fOutD5 = os.path.join(DIR_Annotation, "dump_Mt3.0_annotation.txt");

fInE = os.path.join(DIR_Genome, "Mt3.0_ori.fa");
fOutE = os.path.join(DIR_Genome, "Mt3.0.fa");

attrLst = ["ID", "Name", "Parent", "Note"];
def assemblyOrdered(fIn1, fIn2, fOut):
    fInH1, fInH2, fOutH = open(fIn1, "r"), open(fIn2, "r"), open(fOut, "w");
    for line in fInH1:
        line = line.strip("\n");
        if line != "":
            print >>fOutH, line;
    idLst = [];
    oriDict = {-1:"-", 1:"+"};
    myConv = my.posConv();
    for line in fInH2:
        line = line.strip("\n");
        if line != "":
            ary = line.split("\t");
            if ary[4] != "N":
                chrStr, chrStart, chrStop, ord, type, acc, bacStart, bacStop, \
                    orient, bacLen = ary[0], int(ary[1]), int(ary[2]), ary[3], ary[4], \
                    ary[5], int(ary[6]), int(ary[7]), int(ary[8]), int(ary[9]);
                absStart = myConv.posConversion(1, bacStart, bacStop, chrStart, \
                    chrStop, orient);
                absStop = myConv.posConversion(bacLen, bacStart, bacStop, chrStart, \
                    chrStop, orient);
                if absStart > absStop:
                    absStart, absStop = absStop, absStart;
                assert orient in oriDict;
                chrName, phase, orientConv = \
                    "MtChr"+chrStr, my.phaseBAC[type], oriDict[orient];
                id = acc;
                cnt = 0;
                while id in idLst:
                    cnt += 1;
                    id += "-"+chr(65+cnt);
                idLst.append(id);
                print >>fOutH, "\t".join([chrName, ".", "BAC", str(absStart), \
                    str(absStop), ".", orientConv, ".", \
                    "ID=%s;Name=%s;Note=type(%s)phase(%s)tilingPath(%s:%d..%d)" % \
                    (id, id, type, phase, chrName, chrStart, chrStop)]);
                print >>fOutH, "\t".join([chrName, ".", "tilingBAC", str(chrStart), \
                    str(chrStop), ".", orientConv, ".", \
                    "ID=%s_tiling;Name=%s_tiling;Parent=%s" % (id, id, id)]);
def assemblyUnordered(fIn, fOut):
    fInH, fOutH = open(fIn, "r"), open(fOut, "w");
    cursor = my.conn.cursor(my.MySQLdb.cursors.DictCursor);
    cntO, cntU = dict(), 0;
    for line in fInH:
        line = line.strip("\n");
        if line != "":
            ary = line.split("\t");
            id, length, note = ary[0], int(ary[1]), ary[2];
            qryString = "SELECT chr,acc FROM mt_hapmap.mt3 WHERE acc LIKE '{0}%'";
            cursor.execute(qryString.format(id));
            rst = cursor.fetchall();
            if len(rst) >= 1:
                for row in rst:
                    chr = row["chr"];
                    if chr not in cntO:
                        cntO[chr] = 0;
                    cntO[chr] += 1;
            else:
                cntU += 1;
                print >>fOutH, "\t".join([id, ".", "Unordered_BAC", str(1), \
                    str(length), ".", ".", ".", "ID=%s;Name=%s;Note=%s" % (id,id,note)]);
    fOutH.close();
    print "\n".join(["%d: %d" % (chr, cnt) for chr, cnt in cntO.iteritems()]);
    print "Unordered: %d; Ordered: %d" % (cntU, sum(cntO.values()));
def mergeAssembly(fIn1, fIn2, fOut):
    fInH1, fInH2, fOutH = open(fIn1, "r"), open(fIn2, "r"), open(fOut, "w");
    for line in fInH1:
        line = line.strip("\n");
        if line != "":
            print >>fOutH, line;
    for line in fInH2:
        line = line.strip("\n");
        if line != "":
            print >>fOutH, line;
def checkFormat(fIn, fOut):
    fInH, fOutH = open(fIn, "r"), open(fOut, "w");
    cntLine = 0;
    for line in fInH:
        cntLine += 1;
        line = line.strip("\n");
        if line != "":
            ary = line.split("\t");
            rst = re.match("^chr0(\d)\_pseudomolecule.*", ary[0]);
            if rst != None:
                ary[0] = "MtChr%s" % int(rst.group(1));
            if rst == None:
                rst = re.match("^([A-Za-z]{1,3}[\d]+).*", ary[0]);
                ary[0] = rst.group(1);
            if rst == None:
                print "Abnormal line:" + line;
            attrAry = ary[-1].split(";");
            id, parentId, note = "", "", "";
            for attr in attrAry:
                pair = attr.split("=");
                if len(pair) != 2 or pair[0] not in ["ID", "Parent", "Description"]:
                    print "Invalid line: %d" % cntLine;
                key, value = pair[0], pair[1];
                if key == "ID":
                    id = value;
                elif key == "Parent":
                    parentId = value;
                elif key == "Description":
                    note = value;
            ary[-1] = "ID=%s;Name=%s" % (id,id);
            if parentId != "":
                ary[-1] += ";Parent=%s" % parentId;
            if note != "":
                ary[-1] += ";Note=%s" % note;
            print >>fOutH, "\t".join(ary);
    fOutH.close();
def inferPhase(fIn, fOut1, fOut2):
    fInH, fOutH1 = open(fIn, "r"), open(fOut1, "w");
    mRNAId, orient, lenDict = "", "", dict();
    phaseDict = dict();
    while True:
        line = fInH.readline();
        line = line.strip("\n");
        if not line:
            break;
        else:
            ary = line.split("\t");
            type = ary[2];
            if type == "gene":
                print >>fOutH1, line;
            elif type == "mRNA":
                print >>fOutH1, line;
                attrAry = ary[-1].split(";");
                for attr in attrAry:
                    pair = attr.split("=");
                    assert len(pair) == 2 and pair[0] in attrLst;
                    key, value = pair[0], pair[1];
                    if key == "ID":
                        mRNAId, orient = value, ary[6];
                tmpLineLst = [];
                while True:
                    lineB = fInH.readline();
                    lenLine = len(lineB);
                    lineB = lineB.strip("\n");
                    if lineB != "":
                        aryB = lineB.split("\t");
                        typeB = aryB[2];
                        attrAryB = aryB[-1].split(";");
                        parentId, cdsId = "", "";
                        for attrB in attrAryB:
                            pairB = attrB.split("=");
                            assert len(pairB) == 2 and pairB[0] in attrLst;
                            keyB, valueB = pairB[0], pairB[1];
                            if keyB == "Parent":
                                parentId = valueB;
                            elif keyB == "ID":
                                cdsId = valueB;
                        if parentId == mRNAId:
                            assert typeB in ["CDS", "exon"];
                            tmpLineLst.append(lineB);
                            if typeB == "CDS":
                                #print orient, aryB[6];
                                assert aryB[6] == orient;
                                cdsLen = int(aryB[4]) - int(aryB[3]) + 1;
                                eleAry = cdsId.rsplit("_", 2);
                                if eleAry[0] != mRNAId:
                                    print mRNAId, cdsId;
                                assert eleAry[0]==mRNAId and eleAry[1]=="cds";
                                lenDict[int(aryB[3])] = [int(eleAry[2]), cdsLen];
                        else:
                            mRNAStructure(mRNAId, line, fOutH1, tmpLineLst);
                            fInH.seek(-lenLine, 1);
                            break;
                    else:
                        break;
                lenTot, keysSorted = 0, [];
                if orient == "+":
                    keysSorted = sorted(lenDict.keys());
                else:
                    assert orient == "-";
                    keysSorted = sorted(lenDict.keys(), reverse=True);
                cntCds, cdsOffSet_prev = 0, 0;
                for cdsStart in keysSorted:
                    cntCds += 1;
                    cntCdsI, cdsLen = lenDict[cdsStart][0], lenDict[cdsStart][1];
                    cdsId = "_".join([mRNAId, "cds", str(cntCds)]);
                    cdsIdI = "_".join([mRNAId, "cds", str(cntCdsI)]);
                    if cntCds == 1:
                        cdsOffSet_prev = cntCdsI - cntCds;
                    if cdsOffSet_prev != cntCdsI - cntCds:
                        print cdsIdI+" -> "+cdsId;
                    phaseDict[cdsIdI] = [cdsId, (3 - lenTot % 3) % 3];
                    lenTot += cdsLen;
                if lenTot % 3 != 0:
                    print "\t%s -> potential error" % mRNAId;
                mRNAId, orient, lenDict = "", "", dict();
    with open(fOut2, "w") as fOutH2:
        my.pickle.dump(phaseDict, fOutH2);
    fInH.close();
def mRNAStructure(mRNAId, parentLine, fOutH, tmpLineLst):
    aryParent = parentLine.split("\t");
    cdsMax, cdsMin, orient = int(aryParent[3]), int(aryParent[4]), aryParent[6];
    posCdsLst, posExonLst = [], [];
    featureDict = dict();
    for line in tmpLineLst:
        ary = line.split("\t");
        type, start, stop = ary[2], int(ary[3]), int(ary[4]);
        assert type in ["CDS", "exon"] and orient == ary[6];
        if type == "CDS":
            if start < cdsMin:
                cdsMin = start;
            if stop > cdsMax:
                cdsMax = stop;
            posCdsLst.append([start, stop]);
            featureDict[start] = line;
        else:
            posExonLst.append([start, stop]);
    mygrep = my.grepGene();
    posUtrLst = mygrep.getPosDiff(posExonLst, posCdsLst);
    stopDict = dict();
    for pair in posUtrLst:
        utrStart, utrStop = pair[0], pair[1];
        stopDict[utrStart] = utrStop;
        assert utrStart<=utrStop;
        tag = "";
        if (utrStart>cdsMax and orient=="+") or (utrStop<cdsMin and orient=="-"):
            tag = "three_prime_UTR";
        elif (utrStart>cdsMax and orient=="-") or (utrStop<cdsMin and orient=="+"):
            tag = "five_prime_UTR";
        assert tag != "";
        featureDict[utrStart] = tag;
    cntUtr5, cntUtr3 = 0, 0;
    if orient == "+":
        orderedKeys = sorted(featureDict.keys());
    else:
        orderedKeys = sorted(featureDict.keys(), reverse=True);
    lineOut = "";
    for feStart in orderedKeys:
        tagOrLine = featureDict[feStart];
        if tagOrLine == "five_prime_UTR" or tagOrLine == "three_prime_UTR":
            if tagOrLine == "five_prime_UTR":
                cntUtr5 += 1;
                strTmp, cntUtr = "utr5", cntUtr5;
            else:
                cntUtr3 += 1;
                strTmp, cntUtr = "utr3", cntUtr3;
            lineOut = "\t".join([aryParent[0], aryParent[1], tagOrLine, \
                str(feStart), str(stopDict[feStart]), ".", orient, ".", \
                "ID={0}_{1}_{2};Name={0}_{1}_{2};Parent={0}"
                .format(mRNAId, strTmp, cntUtr)]);
        else:
            lineOut = tagOrLine;
        print >>fOutH, lineOut; 
def addPhase(fIn1, fIn2, fOut):
    fInH1, fOutH = open(fIn1, "r"), open(fOut, "w");
    with open(fIn2, "r") as fInH2:
            phaseDict = my.pickle.load(fInH2);
    while True:
        line = fInH1.readline();
        line = line.strip("\n");
        if not line:
            break;
        else:
            ary = line.split("\t");
            type = ary[2];
            attrAry = ary[-1].split(";");
            if type == "CDS":
                for attr in attrAry:
                    pair = attr.split("=");
                    assert len(pair) == 2 and pair[0] in attrLst;
                    key, value = pair[0], pair[1];
                    if key == "ID":
                        cdsIdI = value;
                        assert cdsIdI in phaseDict;
                        attrAry[attrAry.index(attr)] = \
                            "=".join([key, phaseDict[cdsIdI][0]]);
                        ary[7] = str(phaseDict[cdsIdI][1]);
                ary[-1] = ";".join(attrAry);
        print >>fOutH, "\t".join(ary);
    fInH1.close() and fOutH.close();
def dumpFile(fIn, fOut):
    fOutH = open(fOut, "w");
    with open(fIn, "r") as fInH:
        cnt = 0;
        for line in fInH:
            line = line.strip("\n");
            if line != "":
                ary = line.split("\t");
                assert len(ary) == 9;
                row = { "SeqId":ary[0], "Source":ary[1], "Type":ary[2], "Start":ary[3], "Stop":ary[4], \
                    "Score":ary[5], "Strand":ary[6], "Phase":ary[7], "ID":"", "Parent":"", "Note":""};
                attLst = ary[8].split(";");
                for att in attLst:
                    attPair = att.split("=");
                    assert len(attPair)==2 and attPair[0] in attrLst;
                    row[attPair[0]]=attPair[1];
                del(ary[-1]);
                ary += [row["ID"],row["Parent"],row["Note"]];
                assert len(ary) == 11;
                print >>fOutH, "\t".join(ary);
#LOAD DATA LOCAL INFILE "dump_Mt3.0_annotation.txt" INTO TABLE mt_hapmap.mt3_annotation FIELDS TERMINATED BY "\t"
    
if __name__ == "__main__":  
    #assemblyOrdered(fInA1, fInA2, fOutA);
    #assemblyUnordered(fInB, fOutB);
    #mergeAssembly(fOutA, fOutB, fOutC);
    
    #checkFormat(fInD, fOutD1);
    #inferPhase(fOutD1, fOutD2, fOutD3);
    #addPhase(fOutD2, fOutD3, fOutD4);
    dumpFile(fOutD3, fOutD4);
    #my.rmSeqDesc(fInE, fOutE);

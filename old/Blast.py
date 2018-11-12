import os
import os.path

def getBlastQry(qryPath, mode=1):
    with open(qryPath, 'r') as fQry:
        seqRcdLst = [];
        if mode == 2:
            cntProbeSet = 0;
            cntProbe = 0;
            probeSet_prev = "";
            for seqRcd in SeqIO.parse(fQry, "fasta"):
                tmp = re.split("_", seqRcd.id, 1);
                assert len(tmp) == 2;
                probeSet = tmp[0];
                probe = tmp[1];
                if probeSet == probeSet_prev:
                    cntProbe += 1;
                else:
                    cntProbeSet += 1;
                    cntProbe = 1;
                    probeSet_prev = probeSet;
                if not re.match("\d{3}_\d{3}", seqRcd.id):
                    seqRcd.id = "%03d_%03d" % (cntProbeSet, cntProbe);
                tmpLst = " ".split(seqRcd.description, 1);
                if tmpLst[0] == seqRcd.id:
                    seqRcd.description = tmpLst[1];
                else:
                    seqRcd.description = "%s_%s" % (probeSet, probe);
                seqRcdLst.append(seqRcd);
                #print seqRcd;
        else:
            cnt = 0;
            for seqRcd in SeqIO.parse(fQry, "fasta"):
                cnt += 1;
                if not re.match("\d+", seqRcd.id):
                    seqRcd.id = "%03d" % cnt;
                tmpLst = re.split(" ", seqRcd.description, 1);
                if tmpLst[0] == seqRcd.id:
                    seqRcd.description = tmpLst[1];
                seqRcdLst.append(seqRcd);
    return seqRcdLst;

def blastFilter(fIn, qrySeqRcdLst, fOut, mm_min, mm_max, eValueThresh):
    hitDict = getBlastDict(fIn);
    fOutH = open(fOut,"w");
    for seqRcd in qrySeqRcdLst:
        fOutH.write("%s (%d)\t%s\n" % (seqRcd.id, len(seqRcd), seqRcd.description));
        if not seqRcd.id in hitDict:
            fOutH.write("\tno hit\n");
        else:
            cnt = 0;
            for hsp in hitDict[seqRcd.id]:
                if hsp["expect"] < eValueThresh and len(seqRcd)-hsp["alnLen"]<=mm_max \
                    and len(seqRcd)-hsp["alnLen"]>=mm_min:
                    cnt += 1;
                    fOutH.write("\t{0}\t{1}/{2}\t{3}%\t".format( \
                        "%s:%d..%d"%(hsp["idSbjt"],hsp["startSbjt"],hsp["stopSbjt"]), \
                        hsp["alnLen"], len(seqRcd), hsp["identity"]));
                    fOutH.write("({0} {1} {2} {3} {4} {5})".format( \
                        hsp["startQry"], hsp["stopQry"], hsp["startSbjt"], \
                        hsp["stopSbjt"], hsp["expect"], hsp["bitScore"]));
                    fOutH.write("\n");
            if cnt == 0:
                fOutH.write("\tno hit\n");
    fOutH.close();          

def getBlastDict(fIn):
    hitDict = dict();
    with open(fIn, "r") as fInH:
        for line in fInH:
            eleAry = re.split("\t", line.strip());
            assert len(eleAry)==12;
            eleDict = {"idQry":eleAry[0], "idSbjt":eleAry[1], "identity":float(eleAry[2]), \
                "alnLen":int(eleAry[3])-int(eleAry[4])-int(eleAry[5]), \
                "alnMis":int(eleAry[4]), "alnGap":int(eleAry[5]), "startQry":int(eleAry[6]), "stopQry":int(eleAry[7]), \
                "startSbjt":int(eleAry[8]), "stopSbjt":int(eleAry[9]), "expect":float(eleAry[10]), "bitScore":float(eleAry[11])};
            if not eleDict["idQry"] in hitDict:
                hitDict[eleDict["idQry"]] = [];
            hitDict[eleDict["idQry"]].append(eleDict);
    #print seqLst;
    #print seqMissLst;
    return hitDict;

class blast():
    def __init__(self, qrySeqRcdLst, eValue=5e-2, batchSize=1, \
        optionMask="T", dbOption=1):
        self.blastCmd = "{0}/bin/blastall -p {1} -d {0}/db/{2} -i {3} -o {4} -e {5} -m {6} -a 6 -F '{7}'";
        self.blastDir = os.environ["BLASTDIR"]
        assert self.blastDir != None;
        self.qrySeqRcdLst = qrySeqRcdLst;
        self.eValue = eValue;
        self.batchSize = batchSize;
        self.optionMask = optionMask;
        self.dbOption = dbOption;
        self.seqLst, self.lstLen, self.seqDescLst= [], [], [];
        for seqRcd in self.qrySeqRcdLst:
            self.seqLst.append(seqRcd.id);
            self.seqDescLst.append(seqRcd.description);
            self.lstLen.append(len(seqRcd));
    def blastN(self, fOut):
        blastProgram = "blastn";
        outFormat = 7;
        self.blast(blastProgram, outFormat, fOut);
    def blastP(self, seqRcdLstIn, fOut):
        blastProgram = "tblastn";
        outFormat = 7;
        self.blast(blastProgram, outFormat, fOut);
    def blast(self, blastProgram, outFormat, fOut):
        fOutH = open(fOut, 'w');
        seqRcdLstLst = [];
        seqRcdLst = [];
        cnt = 0;
        for seqRcd in self.qrySeqRcdLst:
            cnt += 1;
            seqRcdLst.append(seqRcd);
            if cnt % self.batchSize == 0:
                seqRcdLstLst.append(seqRcdLst);
                seqRcdLst = [];
        if cnt % self.batchSize != 0:
            seqRcdLstLst.append(seqRcdLst);
        i=0;
        j=0;
        for seqRcdLst in seqRcdLstLst:
            i += 1;
            print "batch %2d: %03d-%03d" % (i, (i-1)*self.batchSize+1, i*self.batchSize);
            ntfQry = tempfile.NamedTemporaryFile(mode='w', delete=False);
            tmpOut = ntfQry.name + ".tmp";
            for seqRcd in seqRcdLst:
                j += 1;
                desc = seqRcd.description;
                if len(desc) > 40: desc=desc[0:40];
                print "\tseq %3d: (%s) %s" % (j,seqRcd.id,desc);
                hdr = seqRcd.id;
                seq = seqRcd.seq.tostring();
                ntfQry.write(">"+hdr);
                ntfQry.write("\n");
                ntfQry.write(seq);
                ntfQry.write("\n");
            ntfQry.close();
            if self.dbOption == 1:
                dbname = "MtChr1-8.fasta";
            elif self.dbOption == 2:
                dbname = "MT3_all_BACs.fas";
            elif self.dbOption == 3:
                dbname = "Mt3.0_cds.fasta";
            elif self.dbOption == 4:
                dbname = "Mt3.0_proteins.fasta";
            else:
                print "error";
                sys.exit(0);
            os.system(self.blastCmd.format(self.blastDir, blastProgram, dbname, \
                ntfQry.name, tmpOut, self.eValue, outFormat, self.optionMask));
            with open(tmpOut, 'r') as fTmpOut:
                for line in fTmpOut:
                    fOutH.write(line);
            #fOut.write(blastCmd.format(blastPath, ntfQry.name, tmpOut) + '\n');
            os.remove(tmpOut);
            os.remove(ntfQry.name);
        fOutH.close();
    def blastExtract(self, fIn, fOut, mm_min, mm_max, evalue):
        seqLst, lstLen, seqDescLst = self.seqLst, self.lstLen, self.seqDescLst;
        seqHitLst, blastDict = [], dict();
        for id in seqLst:
            blastDict[id] = None;
        with open(fIn, "r") as fInH:
            for blastRcd in NCBIXML.parse(fInH):
                qryId = blastRcd.query.split()[0];
                qryIndex = seqLst.index(qryId);
                qryLen = lstLen[qryIndex];
                qryDesc = seqDescLst[qryIndex];
                #print "%s: %s" % (blastRcd.query, qryDesc);
                assert blastRcd.query.split()[0] == seqLst[qryIndex];
                blastRcdDict = dict();
                cntRcdHit = 0;
                for alignment in blastRcd.alignments:
                    sbjctId = re.split(" ", alignment.title, 1)[1];
                    blastAlnDict = dict();
                    cntAlnHit = 0;
                    for hsp in alignment.hsps:
                        if hsp.expect < evalue and len(hsp.query)-hsp.identities<=mm_max \
                                and len(hsp.query)-hsp.identities>=mm_min:
                            hspId = int((hsp.sbjct_start+hsp.sbjct_end)/2);
                            cntRcdHit += 1;
                            cntAlnHit += 1;
                            blastAlnDict[hspId] = { \
                                "query":hsp.query, "sbjct":hsp.sbjct, "match":hsp.match, \
                                "identities":hsp.identities, "positives":hsp.positives, \
                                "query_start":hsp.query_start, \
                                "query_end":hsp.query_end, "sbjct_start":hsp.sbjct_start, \
                                "sbjct_end":hsp.sbjct_end, \
                                "query_length":qryLen, "query_description":qryDesc};
                            #you can do some title matching to extract the subjct id
                    if cntAlnHit > 0:
                        blastRcdDict[sbjctId] = blastAlnDict;
                if cntRcdHit != 0:
                    blastDict[qryId] = blastRcdDict;
        with open(fOut, "w") as fOutH:
            pickle.dump(blastDict, fOutH);
    def blastOutSimple(self, fIn, fOut):
        seqLst, lstLen, seqDescLst = self.seqLst, self.lstLen, self.seqDescLst;
        with open(fIn, "r") as fInH:
            blastDict = pickle.load(fInH);
        fOutH = open(fOut, "w");
        for qryId, blastRcdDict in sorted(blastDict.iteritems()):
            qryIndex = seqLst.index(qryId);
            assert qryId == seqLst[qryIndex];
            qryLen = lstLen[qryIndex];
            qryDesc = seqDescLst[qryIndex];
            fOutH.write("%s\t%s\n" % (qryId, qryDesc));
            for sbjctId, blastAlnDict in sorted(blastRcdDict.iteritems()):
                fOutH.write("\t%s\n" % sbjctId);
                for hspId, hspDict in sorted(blastAlnDict.iteritems()):
                    fOutH.write("\t\t%d..%d\t" % (hspDict["sbjct_start"], hspDict["sbjct_end"])); 
                    fOutH.write("query:%d..%d\t" % (hspDict["query_start"], hspDict["query_end"]));
                    tmp1 = 100 * hspDict["identities"] / len(hspDict["query"]);
                    tmp2 = 100 * hspDict["positives"] / len(hspDict["query"]);
                    fOutH.write("{0}/{1}({2}%)\t{3}/{4}({5}%)\t".format(hspDict["identities"], len(hspDict["query"]), \
                        tmp1, hspDict["positives"], len(hspDict["query"]), tmp2));
                    fOutH.write("\n");
    def blastOut(self, fIn, fOut, fLoc, identyThresh, \
        fSeq, fSeqHit, fSeqNoHit):
        seqLst, lstLen, seqDescLst = self.seqLst, self.lstLen, self.seqDescLst;
        with open(fIn, "r") as fInH:
            blastDict = pickle.load(fInH);
        fOutH = open(fOut, "w");
        fLocH = open(fLoc, "w");
        print >>fLocH, "Gene";
        print >>fLocH, ">DEFL";
        cntHitDict = dict();
        seqEscLst = [];
        for qryId, blastRcdDict in sorted(blastDict.iteritems()):
            qryIndex = seqLst.index(qryId);
            assert qryId == seqLst[qryIndex];
            qryLen = lstLen[qryIndex];
            qryDesc = seqDescLst[qryIndex];
            cntHit = 0;
            if blastRcdDict is None:
                continue;
            for sbjctId, blastAlnDict in sorted(blastRcdDict.iteritems()):
                locRangeDict = dict();
                for loc in blastAlnDict.keys():
                    locRangeDict[loc] = [blastAlnDict[loc]["sbjct_start"], \
                        blastAlnDict[loc]["sbjct_end"]];
                if len(locRangeDict.keys()) > 1000:
                    print "\t%d hits to delete overlaps -> escaped" \
                        % len(locRangeDict.keys());
                    if qryId not in seqEscLst:
                        seqEscLst.append(qryId);
                    continue;
                lstLocAry = [];
                for b in delLocOverlap(locRangeDict):
                    lstLocAry.append([b]);
                if len(lstLocAry) > 300:
                    print "\t%d hits to cluster -> escaped" % len(lstLocAry);
                    if qryId not in seqEscLst:
                        seqEscLst.append(qryId);
                    continue;
                for a in locCluster(lstLocAry):
                    lenHitTot, lenQryTot = 0, 0;
                    qryLocStrLst = [];
                    sbjLocStrLst = [];
                    for locNoOverlap in a:
                        hspDict = blastAlnDict[locNoOverlap];
                        lenHitTot += hspDict["identities"];
                        lenQryTot += len(hspDict["query"]);
                        qryStart, qryEnd, sbjctStart, sbjctEnd = hspDict["query_start"], \
                            hspDict["query_end"], hspDict["sbjct_start"], \
                            hspDict["sbjct_end"];
                        qryLocStrLst.append("%d..%d" % (qryStart, qryEnd));
                        sbjLocStrLst.append("%d..%d" % (sbjctStart, sbjctEnd));
                    ratio = int(lenHitTot/float(qryLen)*100);
                    if ratio > identyThresh:
                        cntHit += 1;
                        desc = qryDesc;
                        if len(desc) > 40: desc=desc[0:40]+"...";
                        if cntHit == 1:
                            cntHitDict[qryId] = 0;
                            print >>fOutH, "%s (%d)\t%s" % (qryId, qryLen, qryDesc);
                            print qryId, desc; 
                        cntHitDict[qryId] += 1;
                        rst = re.findall("(MtChr\d)", sbjctId);
                        assert len(rst) == 1;
                        chrName = rst[0];
                        print >>fOutH, "\t%s\t" % chrName,
                        print >>fOutH, "{0}%\t{1}/{2}/{3}\t".format(\
                            ratio, lenHitTot, lenQryTot, qryLen),
                        print >>fOutH, ";".join(["%s(%s)" % (a,b) \
                            for a,b in zip(qryLocStrLst, sbjLocStrLst)]);
                        print >>fLocH, "%s_%s\t" % (qryId, qryDesc[0:20]) + \
                            ";".join(["%s:%s" % (chrName, j) for j in sbjLocStrLst]) + \
                            "\t"*4;
                        print "\t{0}%".format(ratio) + "\t" + \
                            ";".join(["%s:%s" % (chrName, j) for j in sbjLocStrLst]);
        fOutH.close();
        fLocH.close();
        self.blastSum(cntHitDict, fSeq, fSeqHit, fSeqNoHit);
        for seqEsc in seqEscLst:
            if seqEsc in cntHitDict:
                seqEscLst.pop(seqEscLst.index(seqEsc));
        if len(seqEscLst) > 0:
            print "Also check these escaped queries manually:";
            print " ".join(seqEscLst);
    def blastSum(self, cntHitDict, fSeq, fSeqHit, fSeqNoHit):
        seqRcdHitLst, seqRcdNoHitLst = [], [];
        for seqRcd in self.qrySeqRcdLst:
            if seqRcd.id in cntHitDict:
                seqRcdHitLst.append(seqRcd);
            else:
                seqRcdNoHitLst.append(seqRcd);
        with open(fSeq, "w") as fSeqH:
            SeqIO.write(self.qrySeqRcdLst, fSeqH, "fasta")
        with open(fSeqHit, "w") as fSeqHitH:
            SeqIO.write(seqRcdHitLst, fSeqHitH, "fasta")
        with open(fSeqNoHit, "w") as fSeqNoHitH:
            SeqIO.write(seqRcdNoHitLst, fSeqNoHitH, "fasta")
        print "%d seqs: %d(hits) + %d(no-hits)" % (len(seqRcd), \
            len(seqRcdHitLst), len(seqRcdNoHitLst));
        hitCntSet = set(cntHitDict.values());
        for hitCnt in hitCntSet:
            print "\t%d seqs have %d hits" % (cntHitDict.values().count(hitCnt), \
                hitCnt);        
    def blastParseXML(self, fIn, qrySeqRcdLst, fOut, mm_min, mm_max, evalue):
        fOutH = open(fOut,"w");
        seqHitLst = [];
        seqLst = [];
        lstLen = [];
        seqDescLst = [];
        for seqRcd in qrySeqRcdLst:
            seqLst.append(seqRcd.id);
            seqDescLst.append(seqRcd.description);
            lstLen.append(len(seqRcd));
        cnt = 0;
        with open(fIn, "r") as fInH:
            for blastRcd in NCBIXML.parse(fInH):
                if blastRcd.alignments:
                    pos.append(seqLst.index(blastRcd.query.split()[0]));
                    if blastRcd.query.split()[0] not in seqHitLst:
                        seqHitLst.append(blastRcd.query.split()[0]);
        seqMissSet = set(seqLst) - set(seqHitLst);
        fOutH.write("These sequences have no hits:\n");
        for seqMiss in seqMissSet:
            fOutH.write("%s\t%s\n" % (seqMiss, seqDescLst[seqLst.index(seqMiss)]));
        if len(seqMissSet) == 0:
            fOutH.write("\t(None)\n\n");
        with open(fIn, "r") as fInH:
            for blastRcd in NCBIXML.parse(fInH):
                assert blastRcd.query.split()[0] == seqLst[qryIndex];
                qryId = blastRcd.query.split()[0];
                qryIndex = seqLst.index(qryId);
                qryLen = lstLen[qryIndex];
                qryDesc = seqDescLst[qryIndex];
                print "%s: %s" % (blastRcd.query, qryDesc);
                fOutH.write("%s\t%s\t\n" % (qryId, qryDesc));
                cnt_Hit = 0;
                geneFlag = 0;
                for alignment in blastRcd.alignments:
                    rst = re.match(".*IMGA\|([\w\.]+)\W+(.*)", alignment.title);
                    if rst != None:
                        geneFlag = 1;
                    else:
                        rst = re.match(".*(MtChr\d{1,2}).*", alignment.title);
                    assert rst != None;
                    for hsp in alignment.hsps:
                        if hsp.expect < evalue and len(hsp.query)-hsp.identities<=mm_max \
                                and len(hsp.query)-hsp.identities>=mm_min:
                            cnt_Hit += 1;
                            #if cnt_Hit > 1:
                            fOutH.write("\t");
                            fOutH.write("query:(%d..%d) - %s:%d..%d" % \
                                (hsp.query_start, hsp.query_end, rst.group(1), hsp.sbjct_start, hsp.sbjct_end));
                            fOutH.write("\t{0}/{1}\t{2}/{1}\t".format(hsp.identities, len(hsp.query), hsp.positives));
                            if geneFlag == 1:
                                fOutH.write("\t%s\n" % rst.group(2));
                            #fOutH.write("\t\t%s\n" % (hsp.query));
                            #fOutH.write("\t\t%s\n" % (hsp.match));
                            #fOutH.write("\t\t%s\n" % (hsp.sbjct));
                if cnt_Hit == 0:
                    fOutH.write("\tno hit\n");
        fOutH.close();

suffix = "DEFL";
fQry = os.path.join(my.DIR_In, "in_{0}.fa".format(suffix));
assert os.path.exists(fQry);
DIR_Work = os.path.join(my.DIR_Misc, suffix);
if not os.path.exists(DIR_Work):
    os.mkdir(DIR_Work);
fBlastRst = os.path.join(DIR_Work, "blast.xml");
fBlastDict = os.path.join(DIR_Work, "blast.pck");
fBlastParsed = os.path.join(DIR_Work, "blastParsed.txt");
fOutLocation = os.path.join(DIR_Work, "blastLoc.txt");
fSeq = os.path.join(DIR_Work, "seq.fa");
fSeqHit = os.path.join(DIR_Work, "seqHit_{0}.fa");
fSeqNoHit = os.path.join(DIR_Work, "seqNoHit_{0}.fa");

if __name__ == '__main__':
    option = 1; #set option=2 for probeSet sequence
    seqRcdLst = my.getBlastQry(fQry, option);

    batchSize = 100;
    evalue1 = 5e-5;
    optionMask = 'N';
    dbOption = 2;
    myBlast = my.blast(seqRcdLst, evalue1, batchSize, optionMask, dbOption);
    myBlast.blastN(fBlastRst);
    
    evalue2 = 5e-5;
    mm_min = 0;
    mm_max = 10;
    #myBlast.blastExtract(fBlastRst, fBlastDict, mm_min, mm_max, evalue2);
    identityThresh = 90;
    #myBlast.blastOut(fBlastDict, fBlastParsed, fOutLocation, identityThresh,  fSeq, fSeqHit, fSeqNoHit);


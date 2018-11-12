# -*- coding: utf-8 -*-
import os
import MySQLdb
from Bio import SeqIO
from Bio.Seq import Seq
from BioSQL import BioSeqDatabase

def get_mysql_conn():
    try:
        conn = MySQLdb.connect (host = os.environ["MYSQL_HOST"],
                            user = os.environ["MYSQL_USER"],
                            passwd = os.environ["MYSQL_PW"],
                            db = os.environ["MYSQL_DB"]);
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1]);
        sys.exit(1);
    biosql_server = BioSeqDatabase.open_database( \
        driver="MySQLdb", user=os.environ["MYSQL_USER"], \
        passwd = os.environ["MYSQL_PW"], \
        host = os.environ["MYSQL_HOST"], db=os.environ["MYSQL_DB"]);
    DB = biosql_server["medicago"];

def createSubDB(dbname, desc):
    db = biosql_server.new_database(dbname, description=desc);
    biosql_server.commit();

def loadSeq2DB(seqRcdLst):
    try:
        DB.load(seqRcdLst);
        biosql_server.adaptor.commit();
    except:
        biosql_server.adaptor.rollback();
        raise;

class grepGene():
    def __init__(self, lenMin=400, lenMax=600):
        self.cursor = conn.cursor(MySQLdb.cursors.DictCursor);
        self.lenMin = lenMin;
        self.lenMax = lenMax;
    def getGeneID(self, desc):
        qryString = "SELECT DISTINCT a.ID FROM \
            mt_hapmap.mt3_annotation AS a JOIN mt_hapmap.mt3_annotation AS B \
            ON (a.ID=b.Parent) WHERE ABS(a.Start-a.Stop)>={0} AND \
            ABS(a.Start-a.Stop)<={1} AND b.Note LIKE '%{2}%'";
        self.cursor.execute(qryString.format(self.lenMin, self.lenMax, desc));
        rst = self.cursor.fetchall();
        geneIdLst = [];
        for row in rst:
            if row["ID"] not in geneIdLst:
                geneIdLst.append(row["ID"]);
        return geneIdLst;
    def getPosition(self, geneId, typeLstReq):
        typeLst = ["gene", "mRNA", "exon", "CDS", "intron", "UTR5", "UTR3"];
        for type in typeLstReq:
            assert type in typeLst;
        locLstDict = dict();
        for type in typeLst:
            locLstDict[type] = [];
        qryString = "SELECT SeqID,Start,Stop FROM mt_hapmap.mt3_annotation \
            WHERE ID='{0}'";
        self.cursor.execute(qryString.format(geneId, "gene"));
        rst = self.cursor.fetchall();
        assert len(rst) == 1;
        row = rst[0];
        chrName, geneStart, geneStop = row["SeqID"], row["Start"], row["Stop"];
        loc = "%s:%d..%d" % (chrName, row["Start"], row["Stop"]);
        locLstDict["gene"].append(loc);
        posGeneAry = [[geneStart, geneStop]];
        qryString = "SELECT ID,SeqID,Start,Stop,Strand FROM \
                mt_hapmap.mt3_annotation WHERE \
                ID='{0}.1' AND Type='{1}'";
        self.cursor.execute(qryString.format(geneId, "mRNA", geneStart, \
                geneStop));
        rst = self.cursor.fetchall();
        if len(rst)<1:
            print qryString.format(geneId, "mRNA", geneStart, geneStop);
        assert len(rst) == 1;
        row = rst[0];
        mRNAId = row["ID"];
        orient = row["Strand"];
        loc = "%s:%d..%d" % (chrName, row["Start"], row["Stop"]);
        locLstDict["mRNA"].append(loc);
        posmRNAAry = [[row["Start"], row["Stop"]]];
        qryString = "SELECT SeqID,Start,Stop FROM mt_hapmap.mt3_annotation \
                    WHERE Parent='{0}' AND Type='{1}'";
        self.cursor.execute(qryString.format(mRNAId, "exon"));
        rst = self.cursor.fetchall();
        assert len(rst) >= 1;
        posExonAry = [];
        for row in rst:
            loc = "%s:%d..%d" % (chrName, row["Start"], row["Stop"]);
            locLstDict["exon"].append(loc);
            posExonAry.append([row["Start"], row["Stop"]]);
        posIntronAry = self.getPosDiff(posmRNAAry, posExonAry);
        for posIntron in posIntronAry:
            locLstDict["intron"].append("%s:%d..%d" % \
                (chrName, posIntron[0], posIntron[1]));
        self.cursor.execute(qryString.format(mRNAId, "CDS"));
        rst = self.cursor.fetchall();
        assert len(rst) >= 1;
        posCDSAry = [];
        cdsMin, cdsMax = posmRNAAry[0][1], posmRNAAry[0][0];
        for row in rst:
            loc = "%s:%d..%d" % (chrName, row["Start"], row["Stop"]);
            locLstDict["CDS"].append(loc);
            posCDSAry.append([row["Start"], row["Stop"]]);
            if row["Start"]<cdsMin: cdsMin = row["Start"];
            if row["Stop"]>cdsMax: cdsMax = row["Stop"];
        posUTRAry = self.getPosDiff(posExonAry, posCDSAry);
        for posUTR in posUTRAry:
            assert (posUTR[0]<cdsMin and posUTR[1]<cdsMin) or \
                (posUTR[0]>cdsMax and posUTR[1]>cdsMax);
            if (posUTR[1] < cdsMin and orient == "+") or \
                (posUTR[0] > cdsMax and orient== "-"):
                locLstDict["UTR5"].append("%s:%d..%d" % (chrName, posUTR[0], posUTR[1]));
            elif (posUTR[1] > cdsMax and orient == "+") or \
                (posUTR[0] < cdsMin and orient== "-"):
                locLstDict["UTR3"].append("%s:%d..%d" % (chrName, posUTR[0], posUTR[1]));
        for type in typeLst:
            if type not in typeLstReq:
                del(locLstDict[type]);
        return locLstDict;
    def getPosDiff(self, posAryGlobal, posAryLocal):
        posDiffAry = list(posAryGlobal);
        for posLocal in posAryLocal:
            localStart,localStop = posLocal[0],posLocal[1];
            assert localStart<=localStop;
            for posGlobal in posDiffAry:
                globalStart,globalStop = posGlobal[0],posGlobal[1];
                assert globalStart<=globalStop;
                assert (globalStart<=localStart and globalStop>=localStop) or \
                    (globalStop<localStart) or (globalStart>localStop);
                if globalStart<=localStart and globalStop>=localStop:
                    if globalStart<localStart:
                        posDiffAry.append([globalStart,localStart-1]);
                    if globalStop>localStop:
                        posDiffAry.append([localStop+1,globalStop]);
                    del(posDiffAry[posDiffAry.index(posGlobal)]);
        return posDiffAry;


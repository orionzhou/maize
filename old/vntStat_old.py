# -*- coding: utf-8 -*-
import my;
from numpy import *;

fIn = my.os.path.join(my.DIR_input, "location");
fVnt = my.os.path.join(my.DIR_output, "ret.pck");
fOut = my.os.path.join(my.DIR_stat, 'snpStat.txt');

cutOff = 0.7;
lstAcc = ['HM005_gsnapout8mm', 'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];
	
if __name__ == "__main__":   
	locDict = my.getLoc(fIn);
	fOutH = open(fOut, 'w');
	with open(fVnt, "r") as fVnH:
		rstDict = my.pickle.load(fVnH);
	coverageArySnp, coverageAryIndel = array([]), array([]);
	lstStat = [];
	for loc,rst in sorted(rstDict.iteritems()):
		print "%s\t" % loc;
		locStatDict = dict();
		for key,value in rst.iteritems():
			print "\t%s: (%d)" % (key,len(value.values()));
			coverageDict = dict();
			uniqueAlnDict = dict();
			for pos,row in sorted(value.iteritems()):
				print "\t\t%s(%s/%s)\t" % (pos, row["RefAllele"], row["VarAllele"]);
				for acc in lstAcc:
					accInfo = row[acc];
					print "\t\t\t",
					print "%s\t%s\t%s\t%s\t%.02f\t%s\t%s\t" % (acc, accInfo["Summary"], accInfo["NumReadsWithAllele"], \
						accInfo["Coverage"], float(accInfo["Freq"]), accInfo["UniqAlns"], accInfo["AvgQual"]);
					if acc not in coverageDict:
						coverageDict[acc] = [];
					coverageDict[acc].append(int(accInfo["Coverage"]));
				#print "\t\t\tV:%s" % row["MatchedVariant"];
				#print "\t\t\tR:%s" % row["MatchedReference"];
			#summary about this location
			print "\tSummary";
			for acc in lstAcc:
				if not coverageDict:
					totCov, vntCnt, avgCov= 0, 0, 0;
				else:
					totCov = sum(coverageDict[acc]);
					vntCnt = len(coverageDict[acc]);
					avgCov = totCov / float(vntCnt);
				print "\t\t%s\t%d/%d=%.01f" % (acc, totCov, vntCnt, avgCov);
				oneItem = [loc, acc, key, str(totCov), str(vntCnt), str(avgCov)];
				lstStat.append(oneItem);

	print >>fOutH, "location\tacc\tvariant\ttotalCov\tnumber\tavgCov\t";
	for row in lstStat:
		print >>fOutH, "\t".join(row);
		
		
# -*- coding: utf-8 -*-
import my;
import re;

fBlast = my.os.path.join(my.DIR_output, "blast.out");
fBlastFilterd = my.os.path.join(my.DIR_output, "blast_parsed.out");
fVnt = my.os.path.join(my.DIR_output, "ret.pck");
fOut = my.os.path.join(my.DIR_output, "blast_parsed_with_vnt.out");

cutOff = 0.7;
lstAcc = ['HM005_gsnapout8mm', 'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];
'''lstAcc = ['HM005_gsnapout4mm', 'HM005_gsnapout6mm', 'HM005_gsnapout8mm', \
	'HM029_gsnapout4mm', 'HM029_gsnapout6mm', 'HM029_gsnapout8mm', \
	'HM101_gsnapout4mm', 'HM101_gsnapout6mm', 'HM101_gsnapout8mm'];'''

if __name__ == "__main__":   
	fOutH = open(fOut, "w");
	with open(fVnt, "r") as fVntH:
		rstDict = my.pickle.load(fVntH);
	with open(fBlastFilterd, "r") as fBH:
		for line in fBH:
			rst = re.match("^\t(MtChr\d+:\d+..\d+)\t.*\((.+)\)", line);
			if rst == None:
				fOutH.write(line);
			else:
				line = re.sub("\(.*\)", "", line);
				fOutH.write(line);
				assert rst.group(1) in rstDict;
				strLoc = rst.group(1);
				ary = re.split(" ", rst.group(2));
				assert len(ary) == 6;
				startAln, stopAln, startChr, stopChr, evalue, bitScore = \
					int(ary[0]), int(ary[1]), int(ary[2]), int(ary[3]), float(ary[4]), float(ary[5]);
				flag = 0;
				if startAln > stopAln:
					startAln, stopAln = stopAln, startAln;
					flag += 1;
				if startChr > stopChr:
					startChr, stopChr = stopChr, startChr;
					flag += 1;
				orient = (-1)**flag;
				for key,value in rstDict[strLoc].iteritems():
					vntByAcc = my.snpFilter(value, lstAcc, cutOff, startChr, stopChr, startAln, stopAln, orient);
					for acc, lstVnt in vntByAcc.iteritems():
						cntVnt = 0;
						for vntStr in lstVnt:
							fOutH.write("\t\t");
							if cntVnt == 0:
								fOutH.write("%s (%s)\t" % (acc,key)),
							else:
								fOutH.write("\t"),
							cntVnt += 1;
							fOutH.write("%s\t" % vntStr),
							fOutH.write("\n"),
					#if cntVnt == 0:
						#fOutH.write("\t\tN/A\n");
	fOutH.close();

'''
Created on Aug 3, 2009

@author: rb
'''

import os;
import os.path;
import tempfile;

def main(blastPath, qryPath, outPath, batchSize):
	blastCmd = "{0}/bin/blastall -p blastx -d {0}/db/nr -i {1} -o {2} -e 1e-5 -m 8 -a 6"
	fOut = open(outPath, 'w');
	with open(qryPath, 'r') as fQry:
		hdr = None;
		cnt = 0;
		while hdr != '':
			cnt += 1;
			print("Processing batch %d" % cnt);
			ntfQry = tempfile.NamedTemporaryFile(mode='w', delete=False);
			tmpOut = ntfQry.name + ".out";
			for i in range(batchSize):
				hdr = fQry.readline();
				if (hdr == ''):
					break;
				else:
					seq = fQry.readline();
					ntfQry.write(hdr);
					ntfQry.write(seq);
			ntfQry.close();
			#print(blastCmd.format(blastPath, ntfQry.name, tmpOut));
			os.system(blastCmd.format(blastPath, ntfQry.name, tmpOut));
			with open(tmpOut, 'r') as fTmpOut:
				for line in fTmpOut:
					fOut.write(line);
			#fOut.write(blastCmd.format(blastPath, ntfQry.name, tmpOut) + '\n');
			os.remove(ntfQry.name);
			os.remove(tmpOut);
	fOut.close();
	

if __name__ == '__main__':
	blastPath = "/export/lab/programs/blast";
	qryPath = os.path.join(blastPath, "q/mtgi9not_az");
	outPath = os.path.join(blastPath, "out/mtgi9not_az2.txt");
	batchSize = 20;
	main(blastPath, qryPath, outPath, batchSize);
	
	

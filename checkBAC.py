import my;
import sys;
import os;
import re;

fIn = os.path.join(my.DIR_data, "MT3_all_BACs.fas");
fOut1 = os.path.join(my.DIR_output, "MT3_all_BACs.info");

def getBACinfo(fIn, fOut):
    cnt = 0;
    
  fOutH = open(fOut, "w");
  with open(fIn, "r") as fInH:
    for seqRcd in my.SeqIO.parse(fInH, "fasta"):
      tmp = re.match("%s\W(.*)"%seqRcd.id, seqRcd.description);
      if tmp != None:
        seqRcd.description = tmp.group(1);
      fOutH.write("%s\t%d\t%s\n" % (seqRcd.id,len(seqRcd),seqRcd.description));
      if cnt > 2000000:
        break;
      cnt += 1;
  fOutH.close();
  print "%d BACs processed" % cnt;

def checkBAC(fIn):
  cntA = 0;
  cntB = 0;
  cntC = 0;
  cursor = conn.cursor();
  with open(fIn, "r") as fInH:
    for line in fInH:
      lst = re.split("\t", line.strip());
      print "%s\t%d\t" % (lst[0],int(lst[1])),
      sqlCmd = "SELECT * FROM mt3 WHERE acc LIKE '{0}%'";
      cursor.execute(sqlCmd.format(lst[0]));
      row = cursor.fetchone();
      if row == None:
        print "not found";
        cntB += 1;
      else:
        print "%d\t%d\t%s" % (row[3]-row[2]+1, row[1], row[6]);
        cntA += 1;
      if(cntC > 10000):
        break;
      cntC += 1;
  print "Found:\t%d\nUnfound:\t%d;\nTotal:\t%d" % (cntA, cntB, cntC);
  cursor.close ()
  conn.close();
      
if __name__ == "__main__":   
  #getBACinfo(fIn, fOut1);
  checkBAC(fOut1);

# -*- coding: utf-8 -*-
import my;

fIn = my.os.path.join(my.DIR_input, "location");
fVnt = my.os.path.join(my.DIR_output, "ret.pck");
fOut = my.os.path.join(my.DIR_output, "out.fasta");
   
cutOff = 0.5;
lstAcc = ['HM005_gsnapout8mm', 'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];

if __name__ == "__main__":
    locDict = my.getLoc(fIn);
    with open(fVnt, "r") as fVntH:
        vntDict = my.pickle.load(fVntH);
    mySeqRecover = my.SeqRecover(cutOff, lstAcc, vntDict);
    seqRcdCount = 0;
    seqRcdLst = [];
    for chr, value in sorted(locDict.iteritems()):
        for start, stop in sorted(value.iteritems()):
            seqRcdLst += mySeqRecover.do(chr, start, stop);
    with open(fOut, "w") as fOutH:
        my.SeqIO.write(seqRcdLst, fOutH, "fasta");

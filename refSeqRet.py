# -*- coding: utf-8 -*-
import my;
import sys;
import os.path;

fIn = os.path.join(my.DIR_input, "location");
fOut = os.path.join(my.DIR_output,"refSeq.fasta");

if __name__ == "__main__":
    #retrieve a segment of sequence
    locDict = my.getLoc(fIn);
    myret = my.SeqRet();
    seqRcdLst = [];
    for chr, value in sorted(locDict.iteritems()):
        for start, stop in sorted(value.iteritems()):
            print "MtChr%d:%d..%d\t" % (chr,start,stop),
            seqRcdLst.append(myret.ret(chr, start, stop));
    with open(fOut, "w") as fOutH:
        my.SeqIO.write(seqRcdLst, fOutH, "fasta");

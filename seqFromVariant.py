# -*- coding: utf-8 -*-
import my;

cutOff = 0.3;
fIn = my.os.path.join(my.DIR_input, "location");
fOut = my.os.path.join(my.DIR_output, "out.fasta");
fRet = my.os.path.join(my.DIR_output, "ret.pck");

if __name__ == "__main__":   
    '''lstAcc = ['HM005_gsnapout4mm', 'HM005_gsnapout6mm', 'HM005_gsnapout8mm', \
        'HM029_gsnapout4mm', 'HM029_gsnapout6mm', 'HM029_gsnapout8mm', \
        'HM101_gsnapout4mm', 'HM101_gsnapout6mm', 'HM101_gsnapout8mm'];'''
    lstAcc = ['HM005_gsnapout8mm', 'HM029_gsnapout8mm', 'HM101_gsnapout8mm'];
    seqRecdCount = 0;
    fOutH = open(fOut, "w");
    seqReconstruct = my.SeqReconstruct(cutOff, lstAcc);
    locDict = my.getLoc(fIn, 1);
    #retrieve variants from local file
    variantRet = my.VariantRet();
    #variantRet.ret(locDict, fRet);
    with open(fRet, "r") as fRetH:
        rstDict = my.pickle.load(fRetH);
    for loc,rst in sorted(rstDict.iteritems()):
        print "%s\t\n" % loc,
        for key,value in rst.iteritems():
            print "\t%s: (%d)" % (key,len(value.values())),
            #for pos,row in sorted(value.iteritems()):
            #   print "%01f, " % pos,
            print "\n",
    fOutH.close();
''' seqRecdLst = seqReconstruct.do(lstLoc[0]);
    for seqRecd in seqRecdLst:
        print seqRecd.id + "\t" + seqRecd.description;
    my.SeqIO.write(seqRecdLst, fOutHandle, "fasta");
    if seqRecdCount >= 100000:
        break;
    print "%4d  fasta files have been written" % seqRecdCount;  '''

""" variantQuery = my.VariantQuery();
    variantQuery.open();
    lstAcc = [];
    queryString1 = "SELECT distinct Accession FROM VariantMore";
#   for row in variantQuery.execute(queryString1):
#       if row == None:
#           break;
#      lstAcc = lstAcc + [row[0]];
    variantQuery.close();
    lstAcc = ['HM005_gsnapout8mm', 'HM029_gsnapout8mm', 'HM029_gsnapout6mm'];
    #print lstAcc;
    i = 0;
    fOutHandle = open(fOut, "w");
    seqReconstruct = my.SeqReconstruct_old(cutOff, lstAcc);

    lstString = [];
    with open("location", "r") as fIn:
        for line in fIn:
            pattern = my.re.compile("Chr\d\_\d+\_\d+");
            lstString.append(pattern.match(line.strip()).group());
    #print lstString;
    for str in lstString:
        parts = my.re.split("[\-\_]", str);
        chr, posStart, posStop = int(parts[0][-1]), int(parts[1]), int(parts[2]);         
        seqRecdLst = seqReconstruct.do(chr, posStart, posStop);
        my.SeqIO.write(seqRecdLst, fOutHandle, "fasta");
        i += len(seqRecdLst);
        if i >= 2000:
            break;
    print "%5d  .fasta files have been written" % i;  
    fOutHandle.close();
"""    
"""   
    myCheck = my.ReadsCheck();
    myCheck.open('ReferencesSequences.fasta');
    for recordIn in myCheck.seqIterator:
        parts = my.re.split("[\-\_]",recordIn.id);
        chr, posStart, posStop = int(parts[0][-1]), int(parts[1]), int(parts[2]);       
        seqRecdLst = seqReconstruct.do(chr, posStart, posStop);
        my.SeqIO.write(seqRecdLst, fOutHandle, "fasta");
        i += len(seqRecdLst);
        if i >= 2:
            break;
    print "%5d  .fasta files have been written" % i;  
    fOutHandle.close();
    myCheck.close();
    """

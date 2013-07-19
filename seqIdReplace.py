import my;

fOut = my.os.path.join(my.DIR_input, "RefSeq_rsd.fasta");
if __name__ == "__main__":
    with open(my.DIR_output + "/replace") as Fre:
        lstFind = [];
        lstReplace = [];
        for line in Fre:
            if line.strip() != None:
                parts = line.split("\t");
                lstFind += [parts[0].strip()];
                lstReplace += [parts[1].strip()];
    myCheck = my.ReadsCheck();
    myCheck.open('ReferencesSequences.fasta');
    #print lstFind;
    count = len(lstFind);
    fOutHandle = open(fOut, "w");
    i = 0;
    for recordIn in myCheck.seqIterator:
        for find,replace in zip(lstFind, lstReplace):
            if recordIn.id == find:
                i += 1;
                print find + " -> " + replace;
                recordIn.id = replace;
                recordIn.description = "";
        my.SeqIO.write([recordIn], fOutHandle, "fasta");
    print "\n%d / %d has been replaced" % (i, count);
    myCheck.close();
""" myCheck = my.ReadsCheck();
    myCheck.open('ReferencesSequences_rsd.fasta');
    myCheck.lenCheck();
    myCheck.open('ReferencesSequences_rsd.fasta');
    myCheck.seqCheck();
    myCheck.close();
    myret = my.SeqRet();
"""

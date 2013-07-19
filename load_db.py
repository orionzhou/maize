import my;
import sys;
import os;
import re;

mode = 3;
type = {1:"CDS", 2:"Protein", 3:"BAC_fasta"};
fName = {1:"Mt3.0_cds_20090702_NAMED.fa", \
    2:"Mt3.0_proteins_20090702_NAMED.fa", \
    3:"MT3_all_BACs.fas"};
alphabet = {1:my.IUPAC.unambiguous_dna, 2:my.IUPAC.protein, 3:my.IUPAC.unambiguous_dna};
ctgPattern = "|".join(("chr","AC","CU","fpc"));

def processSeq(mode):
    fIn = os.path.join(my.DIR_data, "Genome_Mt", fName[mode]);
    fInH = open(fIn, "r");
    seqRcdLst = [];
    i = 0;
    for seqRcd in my.SeqIO.parse(fInH, "fasta"):
        if mode == 3:
            tmp = re.match("%s\W(.*)"%seqRcd.id, seqRcd.description);
            if tmp != None:
                seqRcd.description = tmp.group(1);
            print "%s\t%d\t%s" % (seqRcd.id, len(seqRcd), seqRcd.description);
            seqRcdLst.append(seqRcd);
        else:
            seqStr = seqRcd.seq.tostring();
            #seqStr = re.sub('[\W]', '', seqStr);
            lst1 = re.match("^(\w+)\|([\w]+)\.(\d+)$", seqRcd.id);
            lst2 = re.match( "^[\w\|\.]+\W+(.+)\W+((%s)[\w\.]+)\W+(\d+)\-(\d+)\W+(\w)\W+(\w+)\W+(\d+)\W*$" % ctgPattern,  \
                seqRcd.description);
            if lst1==None or lst2==None:
                print "'%s' cannot be identified\n%s" % (seqRcd.id,seqRcd.description);
                sys.exit(1);
            #print "\t".join(lst1.groups());
            #print "\t".join(lst2.groups());
            seqRcdNew = my.SeqRecord(my.Seq(seqStr, alphabet[mode]), id=lst1.group(2)+"."+lst1.group(3));
            seqRcdNew.name = seqRcd.id;
            seqRcdNew.description = lst2.group(1);
            seqRcdNew.annotations['curator'] = lst1.group(1);
            seqRcdNew.annotations['variant'] = lst1.group(3);
            seqRcdNew.annotations['chromosome'] = lst2.group(2);
            seqRcdNew.annotations['start'] = lst2.group(4);
            seqRcdNew.annotations['stop'] = lst2.group(5);
            seqRcdNew.annotations['notes'] = lst2.group(6);
            seqRcdNew.annotations['inst'] = lst2.group(7);
            seqRcdNew.annotations['date'] = lst2.group(8);
            print "%s\t%d\t%s" % (seqRcdNew.id, len(seqRcdNew), seqRcdNew.description);
            #print seqRcdNew;
            seqRcdLst.append(seqRcdNew);
        i += 1;
        if i>10:
            break;
    print "\t %d %s processed" % (i,type[mode]);
    return seqRcdLst;
    
if __name__ == "__main__":
    #my.createSubDB("medicago", "Medicago truncatula");
    seqRcdLst = processSeq(1);
    #load sequences into bio-seqdatabase
    #my.DB.load();
    #my.biosql_server.commit();

import os;
import re;

fo = os.path.join(my.DIR_input, "geneLoc_0102.txt");
retType = ["mRNA", "CDS", "intron", "UTR5", "UTR3"];
geneLst = ["Late nodulin", "Protein kinase", "Peptidase", "Transcriptase", \
    "Leucine-rich repeat", "Oxidase", "Hydrolase", "Helicase", "Transporter"];
geneLst = ["Late nodulin", "Protein kinase", "Peptidase", "Hydrolase", \
    "Oxidase", "Synthase", "Transporter", "Transcriptase"];

lenMax = 3000;
lenMin = 300;
i = 0;
if __name__ == "__main__":   
    mygrep = my.grepGene(lenMin, lenMax);
    fho = open(fo, "w");
    print >>fho, "Gene";
    for geneCatgry in geneLst:
        geneIdLst = mygrep.getGeneID(geneCatgry);
        newId = "%s(%d)" % (re.sub(" ","_",geneCatgry.capitalize()), \
            len(geneIdLst));
        print ">%s" % (newId);
        print >>fho, ">%s" % (newId);
        for geneId in geneIdLst:
            print >>fOutH, "%s\t" % geneId,
            tmpLst = [];
            locLstDict = mygrep.getPosition(geneId, retType);
            for type in retType:
                assert type in locLstDict;
                locLst = locLstDict[type];
                tmpLst.append(";".join(locLst));
            print >>fho, "\t".join(tmpLst);

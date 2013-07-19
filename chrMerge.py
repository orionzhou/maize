import my;
import re;

DIR_chr = my.os.path.join(my.DIR_data, "Genome_Mt");
fOut = my.os.path.join(DIR_chr, "MtChr1-8.fasta");

if __name__ == "__main__":   
    fOutH = open(fOut, "w");
    for i in range(1,9):
        fIn = my.os.path.join(DIR_chr, "MtChr%d"%i);
        with open(fIn, "r") as fInH:
            for line in fInH:
                fOutH.write(line);
    fOutH.close();
    
#formatdb -i MtChr1-8.fasta -p F 

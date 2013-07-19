# -*- coding: utf-8 -*-
import my;
import sys;

fPathIn = my.os.path.join(my.DIR_data, "Variant_Mt", 'variant_table.{0}.10_30_1uniq.txt');
fPathOut1 = my.os.path.join(my.DIR_data, "Variant_Mt", "{0}_snp");
fPathOut2 = my.os.path.join(my.DIR_data, "Variant_Mt", "{0}_indel");
vntChrNameLst = ["MtChloro", "MtMito"] + ["MtChr%d" % i for i in range(1,9)];
vntChrNameLst = [vntChrNameLst[2]];

if __name__ == "__main__":
    #check and split a single variant file into snp and indel file
    print "\t"+"\t".join(["Lines", "snp", "indel"]);
    for vntChrName in vntChrNameLst:
        cntLine = 0;
        fIn = open(fPathIn.format(vntChrName), "r");
        fOut1 = open(fPathOut1.format(vntChrName), "w");
        fOut2 = open(fPathOut2.format(vntChrName), "w");
        count_snp = count_indel = 0;
        pos_prev = 0;
        class_prev = "";
        firstline = fIn.readline();
        fOut1.write(firstline);
        fOut2.write(firstline);
        for line in fIn:
            line = line.strip("\n");
            if line == "":
                break;
            cntLine += 1;
            ary = line.split("\t");
            assert ary[0] == vntChrName;
            if ary[0] != vntChrName:
                ary[0] = vntChrName;
            line = "\t".join(ary);
            class_cur = ary[2];
            if ary[2] == "S":
                count_snp += 1;
                print >>fOut1, line;
                pos_cur = int(ary[3]);
            elif ary[2] == "D" or ary[2] == "I":
                tmp = my.re.match("\D*(\d+)\D*", ary[3]);
                if tmp == None:
                    print "Cannot extract location: %s" % ary[3];
                    sys.exit(1);
                else:
                    count_indel += 1;
                    print >>fOut2, line;
                    pos_cur = int(tmp.group(1));
            else:
                print "Unknow class: %s" % line[0:50];
                sys.exit(1);
            if(class_cur==class_prev and pos_cur<pos_prev):
                print "Not increasing location: %d->%d %s" % (pos_prev, pos_cur, line[0:30]);
                sys.exit(0);
            pos_prev = pos_cur;
            class_prev = class_cur;
        fIn.close();
        print "%s\t%d\t%d\t%d\t" % (vntChrName, cntLine, count_snp, count_indel);
        fOut1.close();
        fOut2.close();
    

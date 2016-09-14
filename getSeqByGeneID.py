import my;
import sys;
import os;
import re;

fIn = my.os.path.join(my.DIR_input, "Germin_Mt_ids.txt");
fOut = my.os.path.join(my.DIR_output, "Germin_Mt.fasta");

if __name__ == "__main__":
    lstId = [];
    lstSeqRecd = []
    with open(fIn, "r") as fInH:
        for line_original in fInH:
            line = line_original.strip();
            tmp = re.match("^([\w\.]+)$", line);
            if tmp != None:
                lstId.append(tmp.group(1));
    for idLocal in lstId:
        seqRecd = my.DB.lookup(primary_id=idLocal);
        seqRecd.description = "";
        lstSeqRecd.append(seqRecd);
        tmp = re.match("^chr(\d+)\_", seqRecd.annotations["chromosome"][0])
        if tmp == None:
            print "cannot determine chromosome: '%s' for %s" % \
                (seqRecd.annotations["chromosome"][0], seqRecd.id);
            sys.exit(1);
        chr = int(tmp.group(1));
        print "%s\tlocation {chr%d:%s-%s}" % (seqRecd.id, chr, \
            seqRecd.annotations["start"][0],seqRecd.annotations["stop"][0]);
    with open(fOut, "w") as fOutH:
        my.SeqIO.write(lstSeqRecd, fOutH, "fasta");

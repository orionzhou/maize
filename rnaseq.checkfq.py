
    fho = open("11.tsv", "w")
    print >>fho, "\t".join([cols[0], "dirf", "read1", "read2"]+cols[1:len(cols)])
    for row in ary:
        rid = row[0]
        print rid
        print >>fhc1, "prefetch %s" % rid
        f_sra = "%s/sra/%s.sra" % (tmp_dir, rid)
        #assert op.isfile(f_sra), "%s not there" % f_sra
        print >>fhc2, "fastq-dump --gzip --split-files -outdir 05.reads %s" % f_sra
        read1, read2 = "%s_1.fastq.gz", "%s_2.fastq.gz"
        f1 = "%s/%s" % (d05, read1)
        f2 = "%s/%s" % (d05, read2)
        #assert op.isfile(f1), "%s not there" % f1
        pe = 0
        if op.isfile(f2):
            pe = 1
        else:
            read2 = ""
        #print >>fho, "\t".join([rid, d05, read1, read2] + row[1:len(row)])

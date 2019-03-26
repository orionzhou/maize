#!/usr/bin/env python
# -*- coding: utf-8 -*-

def get_read_base(aln, qseq, qual, pos):
    for qpos, rpos, rbase in aln:
        if rpos != pos:
            continue
        return [rbase, qseq[qpos], qual[qpos]]
    return ['', '', 0]

def infer_variant(x, minBaseQual):
    aln = x.get_aligned_pairs(matches_only = False, with_seq = True)
    qseq = x.query_sequence
    qual = x.query_qualities
    mms = dict()
    for qpos, rpos, rbase in aln:
        if qpos is not None and rpos is not None and rbase.islower() and qual[qpos] >= minBaseQual:
            mms[qpos] = [rpos+1, 'M', qseq[qpos]]
    #if x.query_name == 'HISEQ15:141:C6VDGANXX:5:1115:10555:79840':
    #    print(mms)
    res = []
    qpos, rpos = 0, x.reference_start
    rbeg = rpos
    vnts = []
    for opt, nbase in x.cigartuples:
        if opt == 0: #M
            for i in range(0, nbase):
                if qpos+i in mms:
                    vnts.append(mms[qpos+i])
            qpos += nbase
            rpos += nbase
        elif opt == 1: #I
            vnts.append([rpos+1, "I", qseq[qpos+1:qpos+nbase+1]])
            qpos += nbase
        elif opt == 2: #D
            vnts.append([rpos+1, "D", nbase])
            rpos += nbase
        elif opt == 3: #N
            res.append([rbeg, rpos, vnts])
            rpos += nbase
            rbeg = rpos
            vnts = []
        elif opt == 4: #S
            qpos += nbase
        elif opt == 7 or opt == 8: #=X
            qpos += 1
            rpos += 1
    if rpos > rbeg:
        res.append([rbeg, rpos, vnts])
    if qpos != x.query_length or rpos != x.reference_end:
        print(qpos, x.query_alignment_end, rpos, x.reference_end, x.cigartuples, aln)
        exit(1)
    return res

def read_variants(fv):
    fhv = open(fv, "r")
    vdic = dict()
    for line in fhv:
        line = line.strip("\n")
        seqid, beg, end, gts = line.split("\t")
        beg, end = int(beg), int(end)
        ref, alt = gts.split(",")
        pos = "%s_%d" % (seqid, int(beg) + 1)
        vdic[pos] = [ref, alt]
    fhv.close()
    return vdic

def bam2bed(args):
    fi, fo = args.bam, args.bed
    min_mapq, min_baseq = args.min_mapq, args.min_baseq
    bam = pysam.AlignmentFile(fi, "rb")
    fho = open(fo, "w")
    for x in bam.fetch():
        if x.is_duplicate or x.is_unmapped or x.is_secondary or x.mapping_quality < min_mapq:
            continue
        sid = x.query_name
        rid, rbeg, rend = x.reference_name, x.reference_start, x.reference_end
        #pair = 'pe' if x.is_paired else 'se'
        #first = 'r1' if x.is_read1 else 'r2'
        res = infer_variant(x, min_baseq)
        #if sid == 'D00635:197:CAC47ANXX:2:2311:18289:97712':
        #    print(res)
        #    exit(1)
        for rbeg, rend, vnts in res:
            vntstr = '.'
            if len(vnts) > 0:
                vntstr = " ".join([":".join(map(str,vnt)) for vnt in vnts])
            fho.write("%s\t%d\t%d\t%s\t%s\n" %
                    (rid, rbeg, rend, sid, vntstr))
        #exit(1)
    fho.close()

def infer_read(rows, fho, fhb):
    rdic, vdic = dict(), dict()
    sid, beg, end, rid = rows[0][0:4]
    locs = set()
    for row in rows:
        assert rid == row[3], "error: " + str(rows)
        if row[4] != '.':
            for rvnt in row[4].split(" "):
                pos, vtype, vnt = rvnt.split(":")
                rdic[pos] = "%s:%s" % (vtype, vnt)
        vpos, vnt, phase = row[7:10]
        assert phase in ['0|1','1|0'], "Unknown phase: %s" % phase
        vdic[vpos] = [vnt, phase]
        loc = "%s:%s" % tuple(row[0:2])
        if loc not in locs:
            locs.add(loc)
            fhb.write("%s\t%s\t%s\t%s\n" % tuple(row[0:4]))
    n0, n1, nunk, nerr = 0, 0, 0, 0
    for pos in vdic:
        if pos in rdic:
            if rdic[pos] == vdic[pos][0]:
                if vdic[pos][1] == '0|1':
                    n1 = n1 + 1
                else:
                    n0 = n0 + 1
            else:
                nunk = nunk + 1
            del rdic[pos]
        else:
            if vdic[pos][1] == '0|1':
                n0 = n0 + 1
            else:
                n1 = n1 + 1
    nerr = len(rdic)
    fho.write("%s\t%d\t%d\t%d\t%d\n" % (rid, n0, n1, nunk, nerr))

def bed_prep(args):
    fi, fo, fb = args.bed, args.out_tsv, args.out_bed
    min_baseq = args.min_baseq
    fhi = open(fi, "r")
    fho = open(fo, "w")
    fhb = open(fb, "w")
    fho.write("rid\tn0\tn1\tnunk\tnerr\n")
    prid = ""
    rows = []
    for line in fhi:
        row = line.strip("\n").split("\t")
        if prid == "":
            rows.append(row)
            prid = row[3]
        elif row[3] == prid:
            rows.append(row)
        else:
            infer_read(rows, fho, fhb)
            rows = [row]
            prid = row[3]
    infer_read(rows, fho, fhb)
    fhi.close()
    fho.close()
    fhb.close()

def bed_summarise(args):
    fr, fb, fo = args.tsv, args.bed, args.out

    fhr = open(fr, "r")
    rdic = dict()
    for line in fhr:
        rid, n0, n1, nunk, nerr = line.strip("\n").split("\t")
        if rid == 'rid':
            continue
        tag = "unk"
        if int(n0) > 0 and int(n1) == 0:
            tag = 'h0'
        elif int(n0) == 0 and int(n1) > 0:
            tag = 'h1'
        elif int(n0) > 0 and int(n1) > 0:
            tag = 'cft'
        rdic[rid] = tag
    fhr.close()

    fhb = open(fb, "r")
    gdic, tdic = dict(), dict()
    for line in fhb:
        row = line.strip("\n").split("\t")
        gid, rid = row[3], row[7]
        if gid not in gdic:
            gdic[gid] = {'n0': 0, 'n1': 0, 'ncft': 0}
        if gid not in tdic:
            tdic[gid] = set()
        if rid not in tdic[gid]:
            tdic[gid].add(rid)
        else:
            continue
        if rid not in rdic:
            print("%s not in read dict" % rid)
        if rdic[rid] == 'h0':
            gdic[gid]['n0'] += 1
        elif rdic[rid] == 'h1':
            gdic[gid]['n1'] += 1
        elif rdic[rid] == 'cft':
            gdic[gid]['ncft'] += 1
    fhb.close()
    
    fho = open(fo, "w")
    fho.write("gid\tn0\tn1\tncft\n")
    for gid in sorted(gdic):
        sdic = gdic[gid]
        n0, n1, ncft = sdic['n0'], sdic['n1'], sdic['ncft']
        if n0 + n1 < 0:
            continue
        fho.write("%s\t%d\t%d\t%d\n" % (gid, n0, n1, ncft))
    fho.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'allele specific expression utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("bam2bed",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'convert bam to bed format using variant coordinates'
    )
    sp1.add_argument('bam', help='input SAM/BAM file')
    sp1.add_argument('bed', help = 'output BED file')
    sp1.add_argument('--min_mapq', default=20, help='min mapping quality')
    sp1.add_argument('--min_baseq', default=20, help='min base quality')
    sp1.add_argument('--vcf', default='/home/springer/zhoux379/data/misc2/mo17vnt/53.vnt.final/61.rna.bed', help='variant file')
    sp1.set_defaults(func = bam2bed)

    sp1 = sp.add_parser("bed_prep",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'infer allele-specific read origin'
    )
    sp1.add_argument('bed_prep', help='input BED (ase) file')
    sp1.add_argument('out_tsv', help = 'output read ASE file (tsv)')
    sp1.add_argument('out_bed', help = 'output read location file (bed)')
    sp1.add_argument('--min_baseq', default=20, help='min base quality')
    sp1.set_defaults(func = bed_prep)

    sp1 = sp.add_parser("bed_summarise",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'summarise allele-specific read counts per gene'
    )
    sp1.add_argument('tsv', help='input tsv (ase) file')
    sp1.add_argument('bed', help = 'input gene-read intersection BED')
    sp1.add_argument('out', help = 'output gene ASE file (tsv)')
    sp1.set_defaults(func = bed_summarise)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


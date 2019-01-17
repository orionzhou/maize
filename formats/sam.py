#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import pysam
import gc

cigar_set = (['M','I','D','N','S','H','P','=','X'])

# SAM columns
sam_QNAME = 0
sam_FLAG = 1
sam_RNAME = 2
sam_POS = 3
sam_MAPQ = 4
sam_CIGAR = 5
sam_RNEXT = 6
sam_PNEXT = 7
sam_TLEN = 8
sam_SEQ = 9
sam_QUAL = 10
sam_TAG = 11

# PSL columns
psl_matches = 0
psl_misMatches = 1
psl_repMatches = 2
psl_nCount = 3
psl_qNumInsert = 4
psl_qBaseInsert = 5
psl_tNumInsert = 6
psl_tBaseInsert = 7
psl_strand = 8
psl_qName = 9
psl_qSize = 10
psl_qStart = 11
psl_qEnd = 12
psl_tName = 13
psl_tSize = 14
psl_tStart = 15
psl_tEnd = 16
psl_blockCount = 17
psl_blockSizes = 18
psl_qStarts = 19
psl_tStarts = 20
psl_seq = 21
psl_empty_line = '0 0 0 0 0 0 0 0 + s 0 0 0 r 0 0 0 0 , , ,'.split()

from maize.apps.base import eprint, sh, mkdir
from maize.formats.base import must_open

def sam2tsv(args):
    sMatch, sMisMatch, sGapOpen, sGapExtend = 2, -3, -5, -2
    sam = pysam.AlignmentFile(args.fi, "r")

    print("qId\tqBeg\tqEnd\tqSrd\tqSize\ttId\ttBeg\ttEnd\ttSrd\ttSize\t" +
            "alnLen\tmatch\tmisMatch\tbaseN\tqNumIns\ttNumIns\tqBaseIns\ttBaseIns\tident\tscore\t" +
            "qLoc\ttLoc")
    for x in sam.fetch():
        if x.is_unmapped:
            continue
        tId, tBeg, tEnd, tSrd, tSize = x.reference_name, x.reference_start, x.reference_end, "+", x.reference_length
        qId, qBeg, qEnd, qSrd, qSize = x.query_name, x.query_alignment_start, x.query_alignment_end, "+", x.query_length
        tBeg += 1
        qBeg += 1
        if x.is_reverse: qSrd = "-"
        if args.paired:
            if x.is_read2:
                qId += ".2"
            else:
                qId += ".1"
        alnLen, match, misMatch, baseN, qLen = 0,0,0,0,0
        qNumIns, tNumIns, qBaseIns, tBaseIns = 0,0,0,0
        for op, nt in x.cigartuples:
            if op == 0 or op == 7 or op == 8: # M=X
                alnLen += nt
                qLen += nt
            elif op == 1: # I
                qNumIns += 1
                qBaseIns += nt
                qLen += nt
            elif op == 4: # S
                qLen += nt
            elif op == 2 or op == 3: # DN
                tNumIns += 1
                tBaseIns += nt
            elif op == 5:
                logging.error("hard clipping: %s -> %s:%d" % (qId, tId, tBeg))
                sys.exit(1)
        if x.has_tag("NM"):
            misMatch = x.get_tag("NM")
        match = alnLen - misMatch
        if qSize == 0:
            qSize = qLen
        #assert qSize == qEnd, "error qSize: %d > %d" % (qSize, qEnd)

        score_match = match * sMatch
        score_misMatch = misMatch * sMisMatch
        numIns = qNumIns + tNumIns
        score_indel = 0
        if numIns >= 1:
            score_indel = sGapOpen + (numIns - 1) * sGapExtend
        score = score_match + score_misMatch + score_indel
        ident = match / (match + misMatch)

        print("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.03f\t%d\t%s\t%s" %
                (qId, qBeg, qEnd, qSrd, qSize,
                 tId, tBeg, tEnd, tSrd, tSize,
                 alnLen, match, misMatch, 0,
                 qNumIns, tNumIns, qBaseIns, tBaseIns, ident, score, '', ''))

def count(args):
    from math import ceil
    winsize, min_qual = args.winsize, args.min_qual
    cdic = dict()
    bam = pysam.AlignmentFile(args.fi, "rb")
    for aln in bam:
        if aln.reference_id == -1:
            seqid = 'unmapped'
            beg = 1
        else:
            seqid = aln.reference_name
            beg = aln.reference_start
        mq = aln.mapping_quality
        if seqid not in cdic:
            cdic[seqid] = dict()
        pbin = ceil(beg / winsize)
        if pbin not in cdic[seqid]:
            cdic[seqid][pbin] = [0, 0]
        if mq >= min_qual:
            cdic[seqid][pbin][0] += 1
        else:
            cdic[seqid][pbin][1] += 1
    bam.close()

    for seqid, idic in sorted(cdic.items()):
        for pbin, cnts in sorted(idic.items()):
            cntstr = "\t".join([str(x) for x in cnts])
            print("%s\t%d\t%s" % (seqid, pbin, cntstr))

def parse_cigar(c,toversion="1.3"):
    # parse CIGAR string
    # TOVERSION: if set to "1.3" the CIGAR is converted forcefully to CIGAR defined
    #          in SAM version 1.3 (i.e. neighbours "X" and "=" are joined and
    #          converted into one region of "M") if CIGAR is given as version 1.4
    #          If it is set to "1.4" no conversion to 1.3 is done!
    r = []
    d = ''
    mismatches_x = 0
    c = c.upper()
    for a in c:
        if a.isdigit():
            d = d + a
        elif a in cigar_set:
            dd = int(d)
            r.append((a,dd))
            if a == 'X':
                mismatches_x = mismatches_x + dd
            d = ''
        else:
            print >>sys.stderr,"ERROR: unknown CIGAR:",c
            sys.exit(1)
    if mismatches_x and toversion == '1.3':
        rr = []
        i = -1
        n = len(r)
        while True:
            i = i + 1
            if i == n:
                r = rr
                break
            elif r[i][0] in ('X', '=', 'M'):
                b = r[i][1]
                for j in xrange(i+1,n):
                    if r[j][0] in ('=','M','X'):
                        b = b + r[j][1]
                    else:
                        i = j - 1
                        break
                rr.append(('M',b))
            else:
                rr.append(r[i])
    return (r,mismatches_x)

def blocks(cigar, ig = 0, use_cigar_13 = True):
    # returns block of matches
    # input is from cigar()
    # NOTE: hard clipping is converted forecfully to soft clipping
    ir = 0 # index on read
    #ig = 0 # index on genome
    rr = [] # on read
    rg = [] # on genome
    match = 0
    mismatch = 0
    mismatch_x = 0
    mismatch_clip = 0
    insert_query = 0
    insert_query_count = 0
    insert_ref = 0
    insert_ref_count = 0
    seq_len = 0 # Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
    #print >>sys.stderr,'cigar:',cigar
    (cig,mismatch_x) = parse_cigar(cigar, toversion = "1.3" if use_cigar_13 else "all")
    mismatch = mismatch_x
    #print >>sys.stderr,"parsed cigar:",cig
    for e in cig:
        if e[0] in ('S','H'):
            ir = ir + e[1] # read
            mismatch = mismatch + e[1]
            mismatch_clip = mismatch_clip + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('I',):
            ir = ir + e[1] # read
            mismatch = mismatch + e[1]
            insert_query = insert_query + e[1]
            insert_query_count = insert_query_count + 1
            seq_len = seq_len + e[1]
        elif e[0] in ('X'):
            ir = ir + e[1] # read
            ig = ig + e[1] # reference/target seq
            mismatch = mismatch + e[1]
            mismatch_x = mismatch_x + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('M','='):
            rr.append((ir,ir+e[1])) # read
            rg.append((ig,ig+e[1])) # reference/target seq
            ir = ir + e[1] # read
            ig = ig + e[1] # reference/target seq
            match = match + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('D','N','P'):
            ig = ig + e[1] # reference/target seq
            insert_ref = insert_ref + e[1]
            insert_ref_count = insert_ref_count + 1
    #print >>sys.stderr,"cigar: query:",rr
    #print >>sys.stderr,"cigar: ref:",rg
    #print >>sys.stderr,"cigar: match:",match
    #print >>sys.stderr,"cigar: mismatch:",mismatch
    #print >>sys.stderr,"cigar: mismatch_clip:",mismatch_clip
    #print >>sys.stderr,"cigar: mismatch_x:",mismatch_x
    return (rr,rg,match,mismatch,mismatch_clip,mismatch_x,insert_ref,insert_ref_count,insert_query,insert_query_count,seq_len)

def get_psl(sam, lens, use_cigar_13=True , replace_string = '', read_sequence=False):
    # USE_CIGAR_13 - If True then the input CIGAR string is in format 1.4 then it will be converted into format 1.3
    #cig, qSize, tSize, tStart, strand):
    # returns PSL coordinates
    # input from blocks()
    #
    #  12. qStart - Alignment start position in query
    #  13. qEnd - Alignment end position in query
    #  18. blockCount - Number of blocks in the alignment (a block contains no gaps)
    #  19. blockSizes - Comma-separated list of sizes of each block
    #  20. qStarts - Comma-separated list of starting positions of each block in query
    #  15. tSize - Target sequence size
    #  16. tStart - Alignment start position in target
    #  17. tEnd - Alignment end position in target
    #  21. tStarts - Comma-separated list of starting positions of each block in target

    psl = None
    if sam and sam[sam_FLAG].isdigit():
        sam_flag = int(sam[sam_FLAG])
        unmapped = True if int(sam[sam_FLAG]) & 0x4  else False
        if (not unmapped) and sam[sam_RNAME] != '*' and sam[sam_CIGAR] != '*' and sam[sam_QNAME] != '*':
            psl = psl_empty_line[:]

            # read sequence length
            psl[psl_tSize] = lens.get(sam[sam_RNAME],0)
            # reference name
            psl[psl_tName] = sam[sam_RNAME]
            # read name
            rname = sam[sam_QNAME].replace(replace_string,'/',1) if replace_string else sam[sam_QNAME]
            psl[psl_qName] = rname
            if (sam_flag & 0x1) and  (sam_flag & 0x40):
                psl[psl_qName] = rname + ".1"
            if (sam_flag & 0x1) and  (sam_flag & 0x80):
                psl[psl_qName] = rname + ".2"

            # strand
            psl[psl_strand] = "-" if int(sam[sam_FLAG]) & 0x10  else '+'

            # start position
            psl[psl_tStart] = int(sam[sam_POS])-1

            (interval_query,interval_ref, match, mismatch, mismatch_clip, mismatch_x,insert_ref,insert_ref_count,insert_query,insert_query_count,seq_len) = blocks(sam[sam_CIGAR], ig = psl[psl_tStart], use_cigar_13 = use_cigar_13)

            # read sequence length
            if sam[sam_SEQ] != '*' and sam[sam_CIGAR].find('H') == -1:
                psl[psl_qSize] = len(sam[sam_SEQ])
            else:
                # • Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
                psl[psl_qSize] = seq_len

            #
            psl[psl_qNumInsert] = insert_query_count
            psl[psl_qBaseInsert] = insert_query
            psl[psl_tNumInsert] = insert_ref_count
            psl[psl_tBaseInsert] = insert_ref

            # extract the mismatches from SAM (using tag NM:i)
            tag_nm_i = [e.partition("NM:i:")[2] for e in sam[sam_TAG:] if e.startswith('NM:i:')] # NM is mismatches per reads
            if not tag_nm_i:
                tag_nm_i = [e.partition("nM:i:")[2] for e in sam[sam_TAG:] if e.startswith('nM:i:')] # nM is not good because but is better than nothing because it is mismatches per fragment and not per read!
            tag_nm_i = int(tag_nm_i[0]) if tag_nm_i else 0
            #if tag_nm_i > float(0.90)*seq_len:
            #    tag_nm_i = 0
            #print >>sys.stderr,"tag NM:i:",tag_nm_i

            # compute the matches and mismatches (include also the clipping as mismatches)
            mis = mismatch_clip + tag_nm_i
            #print >>sys.stderr,"mismatch_clip + tag_nm_i =",mis
            if mis >= mismatch:
                psl[psl_matches] = psl[psl_qSize] - mis
                psl[psl_misMatches] = mis
            else: # probably the tag NM:i are missing???!!!
                psl[psl_matches] = match
                psl[psl_misMatches] = mismatch
            alnsize = sum([e[1]-e[0] for e in interval_query])
            assert alnsize == match, 'alnsize error'
            assert mismatch == mismatch_clip + insert_query, 'q-insert error'
            assert psl[psl_qSize] == match + mismatch_clip + insert_query
            psl[psl_misMatches] = tag_nm_i - insert_query - insert_ref
            psl[psl_matches] = alnsize - psl[psl_misMatches]
            
            #print("------")
            if interval_query:
                #psl[tStart] = boxes[0][1][0]
                psl[psl_tEnd] = interval_ref[-1][1]
                psl[psl_blockCount] = len(interval_query)
                #for i in range(psl[psl_blockCount]):
                #    print("%10d-%10d [%10d]: %10d-%10d [%10d]" % (interval_query[i][0], interval_query[i][1], interval_query[i][1]-interval_query[i][0], interval_ref[i][0], interval_ref[i][1], interval_ref[i][1]-interval_ref[i][0]))

                # this is how is the specification BUT BLAT does not follow the specification!!!
                # NOTE: BLAT _always_ gives the coordinates as everything is mapped on the forwward strand
                # even that it is mapped on the reverse strand
                #
                if psl[psl_strand] == "+":
                    psl[psl_qStart] = interval_query[0][0]
                    psl[psl_qEnd] = interval_query[-1][1]
                #    psl[psl_blockSizes] = ','.join([str(e[1]-e[0]) for e in interval_query])+','
                #    psl[psl_qStarts] = ','.join([str(e[0]) for e in interval_query])+','
                #    psl[psl_tStarts] = ','.join([str(e[0]) for e in interval_ref])+','
                elif psl[psl_strand] == "-":
                    psl[psl_qStart] = psl[psl_qSize] - interval_query[-1][1]
                    psl[psl_qEnd] = psl[psl_qSize] - interval_query[0][0]
                #    psl[psl_blockSizes] = ','.join([str(e[1]-e[0]) for e in interval_query[::-1]])+','
                #    psl[psl_qStarts] = ','.join([str(psl[psl_qSize]-e[1]) for e in interval_query[::-1]])+','
                #    psl[psl_tStarts] = ','.join([str(psl[psl_tSize]-e[0]-1) for e in interval_ref[::-1]])+','

                psl[psl_blockSizes] = ','.join([str(e[1]-e[0]) for e in interval_query])+','
                psl[psl_qStarts] = ','.join([str(e[0]) for e in interval_query])+','
                psl[psl_tStarts] = ','.join([str(e[0]) for e in interval_ref])+','

            if read_sequence:
                if len(psl) < psl_seq + 1:
                    psl.append(sam[sam_SEQ])

            psl = map(str,psl)

    return psl

def getlines(a_filename):
    # it gives chunks
    fin = None
    if a_filename == '-':
        fin = sys.stdin
    else:
        fin = open(a_filename,'r')
    header = dict()
    first = True
    while True:
        lines = fin.readlines(10**8)
        if not lines:
            break
        gc.disable()
        lines = [line.rstrip('\r\n').split('\t') for line in lines if line.rstrip('\r\n')]
        gc.enable()
        for line in lines:
            if line[0].startswith('@'):
                if line[0].startswith('@SQ') and line[1].startswith('SN:') and line[2].startswith('LN:'):
                    k = line[1][3:]
                    v = int(line[2][3:])
                    header[k] = v
                else:
                    pass
            else:
                if first:
                    first = False
                    yield header
                    header = None
                yield line
    if first and header:
        yield header
    fin.close()

def sam2psl(args):
    file_in, file_ou = args.sam, args.psl
    use_cigar_13 = not args.skip_conversion_cigar_13
    read_sequence = args.read_sequence
    replace_string = args.replace_reads_ids if args.replace_reads_ids else ''

    fou = None
    if file_ou == '-':
        fou = sys.stdout
    else:
        fou = open(file_ou,'w')

    # PSL data
    psl = []
    psl_empty_line = ['0']*21
    # processing
    i = 0
    size_lines = 10**6
    lengths = None
    for line in getlines(file_in):
        if i == 0:
            lengths = line
            i = i + 1
            continue
        i = i + 1
        temp = get_psl(line, lengths, use_cigar_13, replace_string, read_sequence)
        # saving
        if temp:
            psl.append('\t'.join(temp)+'\n')
            #print >>sys.stderr, '\t'.join(line)
            #print >>sys.stderr, '\t'.join(temp)
            #print >>sys.stderr, "-----------------------------------------------"
            if i > size_lines:
                fou.writelines(psl)
                psl = []
    if psl:
        fou.writelines(psl)
    fou.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'sam utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("2tsv", help = "sam -> tsv")
    sp1.add_argument('fi', help = 'input *.sam or *.bam file')
    sp1.add_argument('--paired', action = "store_true", help = 'paired end input ')
    sp1.set_defaults(func = sam2tsv)

    sp1 = sp.add_parser("2psl", help = "sam -> psl",
                        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument("sam", help="The input file in SAM format.")
    sp1.add_argument("psl", help="The output file in PSL format.")
    sp1.add_argument("--skip-conversion-cigar-1.3","-4", action = "store_true",
                     dest = "skip_conversion_cigar_13",
                     help="By default if the CIGAR strings in the input SAM file are in the format defined in "+
                     "SAM version 1.4 (i.e. there are 'X' and '=') then the CIGAR string will be "+
                     "first converted into CIGAR string, which is described in SAM version 1.3, "+
                     "(i.e. there are no 'X' and '=' which are replaced with 'M') and afterwards into PSL format.")
    sp1.add_argument("--read-seq","-s", action = "store_true", dest = "read_sequence",
                      help = "It adds to the PSL output as column 22, the sequence of the read. "+
                     "This is not anymore a valid PSL format.")
    sp1.add_argument("--replace-read-ids","-r", dest = "replace_reads_ids",
                     help = "In the reads ids (also known as query name in PSL) "+
                     "the string specified here will be replaced with '/' "+
                     "(which is used in Solexa for /1 and /2).")
    sp1.set_defaults(func = sam2psl)

    sp1 = sp.add_parser("count",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'count reads in bins'
    )
    sp1.add_argument('fi', help='input *.sam or *.bam file')
    sp1.add_argument('--winsize', default=100000, choices=[100000, 1000000], help='window size')
    sp1.add_argument('--min_qual', default=20, choices=[0, 10, 20], help='min read mapping quality')
    sp1.set_defaults(func = count)
 
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()



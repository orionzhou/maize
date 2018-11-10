#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging
import string

from collections import defaultdict
from itertools import product, combinations

from maize.formats.blast import BlastLine
from maize.formats.fasta import Fasta
from maize.formats.bed import Bed
from maize.formats.base import must_open, BaseFile
from maize.utils.grouper import Grouper
from maize.utils.cbook import gene_name
from maize.compara.synteny import AnchorFile, check_beds
from maize.apps.base import glob, need_update, sh, mkdir

class OMGFile (BaseFile):

    def __init__(self, filename):
        super(OMGFile, self).__init__(filename)
        fp = open(filename)
        inblock = False
        components = []
        component = []
        for row in fp:
            if inblock:
                atoms = row.split()
                natoms = len(atoms)
                assert natoms in (0, 7)
                if natoms:
                    gene, taxa = atoms[0], atoms[5]
                    component.append((gene, taxa))
                else:
                    inblock = False
                    components.append(tuple(component))

            if row.strip().startswith("---"):
                inblock = True
                component = []

        if inblock:
            components.append(tuple(component))
        self.components = components

    def best(self):
        bb = set()
        for component in self.components:
            size = len(component)
            if size > 1:
                bb.add(component)
        return bb

def get_weights(weightsfiles=None):
    if weightsfiles is None:
        weightsfiles = glob("*.weights")

    weights = defaultdict(list)
    for row in must_open(weightsfiles):
        a, b, c = row.split()
        weights[a].append((a, b, c))
    return weights

def get_edges(weightsfiles=None):
    if weightsfiles is None:
        weightsfiles = glob("*.weights")

    edges = {}
    for row in must_open(weightsfiles):
        a, b, c = row.split()
        c = int(c)
        edges[(a, b)] = c
        edges[(b, a)] = c
    return edges

def get_info():
    infofiles = glob("*.info")
    info = {}
    for row in must_open(infofiles):
        a = row.split()[0]
        info[a] = row.rstrip()
    return info

def pairwise_distance(a, b, threadorder):
    d = 0
    for x, y in zip(a, b)[:-1]:  # Last column not used
        x, y = x.strip("|"), y.strip("|")
        if "." in (x, y):
            dd = 50
        else:
            xi, x = threadorder[x]
            yi, y = threadorder[y]
            dd = min(abs(xi - yi), 50)
        d += dd
    return d

def insert_into_threaded(atoms, threaded, threadorder):
    min_idx, min_d = 0, 1000
    for i, t in enumerate(threaded):
        # calculate distance
        d = pairwise_distance(atoms, t, threadorder)
        if d < min_d:
            min_idx = i
            min_d = d

    i = min_idx
    t = threaded[i]
    threaded.insert(i, atoms)
    logging.debug("Insert {0} before {1} (d={2})".format(atoms, t, min_d))

def sort_layout(thread, listfile, column=0):
    """
    Sort the syntelog table according to chromomomal positions. First orient the
    contents against threadbed, then for contents not in threadbed, insert to
    the nearest neighbor.
    """
    from maize.formats.base import DictFile

    outfile = listfile.rsplit(".", 1)[0] + ".sorted.list"
    threadorder = thread.order
    fw = open(outfile, "w")
    lt = DictFile(listfile, keypos=column, valuepos=None)
    threaded = []
    imported = set()
    for t in thread:
        accn = t.accn
        if accn not in lt:
            continue

        imported.add(accn)
        atoms = lt[accn]
        threaded.append(atoms)

    assert len(threaded) == len(imported)

    total = sum(1 for x in open(listfile))
    logging.debug("Total: {0}, currently threaded: {1}".format(total, len(threaded)))
    fp = open(listfile)
    for row in fp:
        atoms = row.split()
        accn = atoms[0]
        if accn in imported:
            continue
        insert_into_threaded(atoms, threaded, threadorder)

    for atoms in threaded:
        print >> fw, "\t".join(atoms)

    fw.close()
    logging.debug("File `{0}` sorted to `{1}`.".format(outfile, thread.filename))

def layout(args):
    """
    %prog layout omgfile taxa

    Build column formatted gene lists after omgparse(). Use species list
    separated by comma in place of taxa, e.g. "BR,BO,AN,CN"
    """

    omgfile, taxa = args.omgfile, args.taxa
    listfile = omgfile.rsplit(".", 1)[0] + ".list"
    taxa = taxa.split(",")
    ntaxa = len(taxa)
    fw = open(listfile, "w")

    data = []
    fp = open(omgfile)
    for row in fp:
        genes, idxs = row.split()
        row = ["."] * ntaxa
        genes = genes.split(",")
        ixs = [int(x) for x in idxs.split(",")]
        for gene, idx in zip(genes, ixs):
            row[idx] = gene
        txs = ",".join(taxa[x] for x in ixs)
        print >> fw, "\t".join(("\t".join(row), txs))
        data.append(row)

    coldata = zip(*data)
    ngenes = []
    for i, tx in enumerate(taxa):
        genes = [x for x in coldata[i] if x != '.']
        genes = set(x.strip("|") for x in genes)
        ngenes.append((len(genes), tx))

    details = ", ".join("{0} {1}".format(a, b) for a, b in ngenes)
    total = sum(a for a, b in ngenes)
    s = "A list of {0} orthologous families that collectively".format(len(data))
    s += " contain a total of {0} genes ({1})".format(total, details)
    print >> sys.stderr, s

    fw.close()
    lastcolumn = ntaxa + 1
    cmd = "sort -k{0},{0} {1} -o {1}".format(lastcolumn, listfile)
    sh(cmd)

    logging.debug("List file written to `{0}`.".format(listfile))
    sort = args.sort
    if sort:
        thread = Bed(sort)
        sort_layout(thread, listfile)

def omgparse(args):
    """
    %prog omgparse work

    Parse the OMG outputs to get gene lists.
    """
    work = args.work
    omgfiles = glob(op.join(work, "gf*.out"))
    for omgfile in omgfiles:
        omg = OMGFile(omgfile)
        best = omg.best()
        for bb in best:
            genes, taxa = zip(*bb)
            print( "\t".join((",".join(genes), ",".join(taxa))) )

def group(args):
    """
    %prog group anchorfiles

    Group the anchors into ortho-groups. Can input multiple anchor files.
    """

    anchorfiles = args.anchorfiles
    groups = Grouper()
    for anchorfile in anchorfiles:
        ac = AnchorFile(anchorfile)
        for a, b, idx in ac.iter_pairs():
            groups.join(a, b)

    logging.debug("Created {0} groups with {1} members.".\
                  format(len(groups), groups.num_members))

    outfile = args.outfile
    fw = must_open(outfile, "w")
    for g in groups:
        print >> fw, ",".join(sorted(g))
    fw.close()

    return outfile

def omg(args):
    """
    %prog omg weightsfile

    Run Sankoff's OMG algorithm to get orthologs. Download OMG code at:
    <http://137.122.149.195/IsbraSoftware/OMGMec.html>

    This script only writes the partitions, but not launch OMGMec. You may need to:

    $ parallel "java -cp ~/code/OMGMec TestOMGMec {} 4 > {}.out" ::: work/gf?????

    Then followed by omgparse() to get the gene lists.
    """

    weightsfiles = args.weightsfiles
    groupfile = group(weightsfiles + ["--outfile=groups"])

    weights = get_weights(weightsfiles)
    info = get_info()

    fp = open(groupfile)

    work = "work"
    mkdir(work)
    for i, row in enumerate(fp):
        gf = op.join(work, "gf{0:05d}".format(i))
        genes = row.rstrip().split(",")

        fw = open(gf, "w")
        contents = ""
        npairs = 0
        for gene in genes:
            gene_pairs = weights[gene]
            for a, b, c in gene_pairs:
                if b not in genes:
                    continue

                contents += "weight {0}".format(c) + '\n'
                contents += info[a] + '\n'
                contents += info[b] + '\n\n'
                npairs += 1

        header = "a group of genes  :length ={0}".format(npairs)
        print >> fw, header
        print >> fw, contents

        fw.close()

def geneinfo(bed, order, genomeidx, ploidy):
    bedfile = bed.filename
    p = bedfile.split(".")[0]
    idx = genomeidx[p]
    pd = ploidy[p]
    infofile = p + ".info"

    if not need_update(bedfile, infofile):
        return infofile

    fwinfo = open(infofile, "w")

    for s in bed:
        chr = "".join(x for x in s.seqid if x in string.digits)
        try:
            chr = int(chr)
        except ValueError:
            chr = "0"

        print >> fwinfo, "\t".join(str(x) for x in \
                    (s.accn, chr, s.start, s.end, s.strand, idx, pd))
    fwinfo.close()

    logging.debug("Update info file `{0}`.".format(infofile))

    return infofile

def omgprepare(args):
    """
    %prog omgprepare ploidy anchorsfile blastfile

    Prepare to run Sankoff's OMG algorithm to get orthologs.
    """
    from maize.formats.blast import cscore
    from maize.formats.base import DictFile

    ploidy, anchorfile, blastfile = args.ploidy, args.anchorfile, args.blastfile
    norbh = args.norbh
    pctid = args.pctid
    cs = args.cscore
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    fp = open(ploidy)
    genomeidx = dict((x.split()[0], i) for i, x in enumerate(fp))
    fp.close()

    ploidy = DictFile(ploidy)

    geneinfo(qbed, qorder, genomeidx, ploidy)
    geneinfo(sbed, sorder, genomeidx, ploidy)

    pf = blastfile.rsplit(".", 1)[0]
    cscorefile = pf + ".cscore"
    cscore([blastfile, "-o", cscorefile, "--cutoff=0", "--pct"])
    ac = AnchorFile(anchorfile)
    pairs = set((a, b) for a, b, i in ac.iter_pairs())
    logging.debug("Imported {0} pairs from `{1}`.".format(len(pairs), anchorfile))

    weightsfile = pf + ".weights"
    fp = open(cscorefile)
    fw = open(weightsfile, "w")
    npairs = 0
    for row in fp:
        a, b, c, pct = row.split()
        c, pct = float(c), float(pct)
        c = int(c * 100)
        if (a, b) not in pairs:
            if norbh:
                continue
            if c < cs:
                continue
            if pct < pctid:
                continue
            c /= 10  # This severely penalizes RBH against synteny

        print >> fw, "\t".join((a, b, str(c)))
        npairs += 1
    fw.close()

    logging.debug("Write {0} pairs to `{1}`.".format(npairs, weightsfile))

def make_ortholog(blocksfile, rbhfile, orthofile):
    from maize.formats.base import DictFile

    # Generate mapping both ways
    adict = DictFile(rbhfile)
    bdict = DictFile(rbhfile, keypos=1, valuepos=0)
    adict.update(bdict)

    fp = open(blocksfile)
    fw = open(orthofile, "w")
    nrecruited = 0
    for row in fp:
        a, b = row.split()
        if b == '.':
            if a in adict:
                b = adict[a]
                nrecruited += 1
                b += "'"
        print >> fw, "\t".join((a, b))

    logging.debug("Recruited {0} pairs from RBH.".format(nrecruited))
    fp.close()
    fw.close()

def ortholog(args):
    """
    %prog ortholog species_a species_b

    Run a sensitive pipeline to find orthologs between two species a and b.
    The pipeline runs LAST and generate .lifted.anchors.

    `--full` mode would assume 1-to-1 quota synteny blocks as the backbone of
    such predictions. Extra orthologs will be recruited from reciprocal best
    match (RBH).
    """
    from maize.apps.align import last as last_main
    from maize.compara.blastfilter import main as blastfilter_main
    from maize.compara.quota import main as quota_main
    from maize.compara.synteny import scan, mcscan, liftover
    from maize.formats.blast import cscore, filter

    a, b = args.species_a, args.species_b
    dbtype = args.dbtype
    suffix = ".cds" if dbtype == "nucl" else ".pep"
    abed, afasta = a + ".bed", a + suffix
    bbed, bfasta = b + ".bed", b + suffix
    ccscore = args.cscore
    quota = args.quota
    dist = "--dist={0}".format(args.dist)

    aprefix = afasta.split(".")[0]
    bprefix = bfasta.split(".")[0]
    pprefix = ".".join((aprefix, bprefix))
    qprefix = ".".join((bprefix, aprefix))
    last = pprefix + ".last"
    if need_update((afasta, bfasta), last):
        last_main([bfasta, afasta], dbtype)

    if a == b:
        lastself = last + ".P98L0.inverse"
        if need_update(last, lastself):
            filter([last, "--hitlen=0", "--pctid=98", "--inverse", "--noself"])
        last = lastself

    filtered_last = last + ".filtered"
    if need_update(last, filtered_last):
        if args.no_strip_names:
            blastfilter_main([last, "--cscore={0}".format(ccscore), "--no_strip_names"])
        else:
            blastfilter_main([last, "--cscore={0}".format(ccscore)])

    anchors = pprefix + ".anchors"
    lifted_anchors = pprefix + ".lifted.anchors"
    pdf = pprefix + ".pdf"
    if not args.full:
        if need_update(filtered_last, lifted_anchors):
            if args.no_strip_names:
                scan([filtered_last, anchors, dist,
                        "--liftover={0}".format(last), "--no_strip_names"])
            else:
                scan([filtered_last, anchors, dist,
                        "--liftover={0}".format(last)])
        if quota:
            quota_main([lifted_anchors,
                        "--quota={0}".format(quota), "--screen"])
        if need_update(anchors, pdf):
            from maize.graphics.dotplot import dotplot_main
            dargs = [anchors]
            if args.nostdpf:
                dargs += ["--nostdpf", "--skipempty"]
            dotplot_main(dargs)
        return

    if need_update(filtered_last, anchors):
        if args.no_strip_names:
            scan([filtered_last, anchors, dist, "--no_strip_names"])
        else:
            scan([filtered_last, anchors, dist])

    ooanchors = pprefix + ".1x1.anchors"
    if need_update(anchors, ooanchors):
        quota_main([anchors, "--quota=1:1", "--screen"])

    lifted_anchors = pprefix + ".1x1.lifted.anchors"
    if need_update((last, ooanchors), lifted_anchors):
        if args.no_strip_names:
            liftover([last, ooanchors, dist, "--no_strip_names"])
        else:
            liftover([last, ooanchors, dist])

    pblocks = pprefix + ".1x1.blocks"
    qblocks = qprefix + ".1x1.blocks"
    if need_update(lifted_anchors, [pblocks, qblocks]):
        mcscan([abed, lifted_anchors, "--iter=1", "-o", pblocks])
        mcscan([bbed, lifted_anchors, "--iter=1", "-o", qblocks])

    rbh = pprefix + ".rbh"
    if need_update(last, rbh):
        cscore([last, "-o", rbh])

    portho = pprefix + ".ortholog"
    qortho = qprefix + ".ortholog"
    if need_update([pblocks, qblocks, rbh], [portho, qortho]):
        make_ortholog(pblocks, rbh, portho)
        make_ortholog(qblocks, rbh, qortho)

def tandem_main(blast_file, cds_file, bed_file, N=3, P=50, is_self=True, \
    evalue=.01, strip_name=".", ofile=sys.stderr, genefam=False):

    if genefam:
        N = 1e5

    # get the sizes for the CDS first
    f = Fasta(cds_file)
    sizes = dict(f.itersizes())

    # retrieve the locations
    bed = Bed(bed_file)
    order = bed.order

    if is_self:
        # filter the blast file
        g = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len)*P/100.:
                continue

            query = gene_name(b.query, strip_name)
            subject = gene_name(b.subject, strip_name)
            qi, q = order[query]
            si, s = order[subject]

            if abs(qi - si) <= N and b.evalue <= evalue:
                if genefam:
                    g.join(query, subject)
                elif q.seqid == s.seqid:
                    g.join(query, subject)

    else:
        homologs = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len)*P/100.:
                continue
            if b.evalue > evalue:
                continue

            query = gene_name(b.query, strip_name)
            subject = gene_name(b.subject, strip_name)
            homologs.join(query, subject)

        if genefam:
            g = homologs
        else:
            g = Grouper()
            for i, atom in enumerate(bed):
                for x in range(1, N+1):
                    if all([i-x >= 0, bed[i-x].seqid == atom.seqid, \
                        homologs.joined(bed[i-x].accn, atom.accn)]):
                        leni = sizes[bed[i].accn]
                        lenx = sizes[bed[i-x].accn]
                        if abs(leni - lenx) > max(leni, lenx)*(1-P/100.):
                            continue
                        g.join(bed[i-x].accn, atom.accn)

    # dump the grouper
    fw = must_open(ofile, "w")
    ngenes, nfamilies = 0, 0
    families = []
    for group in sorted(g):
        if len(group) >= 2:
            print >>fw, ",".join(sorted(group))
            ngenes += len(group)
            nfamilies += 1
            families.append(sorted(group))

    longest_family = max(families, key=lambda x: len(x))

    # generate reports
    print >>sys.stderr, "Proximal paralogues (dist=%d):" % N
    print >>sys.stderr, "Total %d genes in %d families" % (ngenes, nfamilies)
    print >>sys.stderr, "Longest families (%d): %s" % (len(longest_family),
        ",".join(longest_family))

    return families

def tandem(args):
    """
    %prog tandem blast_file cds_file bed_file [options]

    Find tandem gene clusters that are separated by N genes, based on filtered
    blast_file by enforcing alignments between any two genes at least 50%
    (or user specified value) of either gene.

    pep_file can also be used in same manner.
    """
    blast_file, cds_file, bed_file = args.blast_file, args.cds_file, args.bed_file
    N = args.tandem_Nmax
    P = args.percent_overlap
    is_self = not args.not_self
    sep = args.sep
    ofile = args.out_file

    tandem_main(blast_file, cds_file, bed_file, N=N, P=P, is_self=is_self, \
        evalue=args.evalue, strip_name=sep, ofile=ofile, genefam=args.genefam)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'catalog utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("tandem", 
            help = 'identify tandem gene groups within certain distance',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('blast_file', help = 'blast file')
    sp1.add_argument('cds_file', help = 'cds file')
    sp1.add_argument('bed_file', help = 'bed file')
    sp1.add_argument('out_file', help = 'out file')
    sp1.add_argument("--tandem_Nmax", dest="tandem_Nmax", type=int, default=3,
               help="merge tandem genes within distance")
    sp1.add_argument("--percent_overlap", type=int, default=50,
               help="tandem genes have >=x%% aligned sequence, x=0-100")
    sp1.add_argument("--evalue", type=float, default=0.01, help = 'evalue')
    sp1.add_argument("--not_self", action="store_true",
                 help="provided is not self blast file")
    sp1.add_argument("--strip_gene_name", dest="sep", default=".",
               help="strip alternative splicing. Use None for no stripping")
    sp1.add_argument("--genefamily", dest="genefam", action="store_true",
                 help="compile gene families based on similarity")
    sp1.set_defaults(func = tandem)

    sp1 = sp.add_parser("ortholog", 
            help = 'run a combined synteny and RBH pipeline to call orthologs',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('species_a', help = 'species A')
    sp1.add_argument('species_b', help = 'species B')
    sp1.add_argument("--dbtype", default="nucl", choices=("nucl", "prot"),
                 help="Molecule type of subject database")
    sp1.add_argument("--full", action="store_true",
                 help="Run in full mode, including blocks and RBH")
    sp1.add_argument("--cscore", default=0.7, type=float,
                 help="C-score cutoff")
    sp1.add_argument("--dist", default=20, type=int,
                 help="Extent of flanking regions to search")
    sp1.add_argument("--quota", help="Quota align parameter")
    sp1.add_argument("--nostdpf", action="store_true",
            help="Do not standardize contig names")
    sp1.add_argument("--no_strip_names", action="store_true",
            help="Do not strip alternative splicing "
            "(e.g. At5g06540.1 -> At5g06540)")
    sp1.set_defaults(func = ortholog)

    sp1 = sp.add_parser("group", 
            help = 'cluster the anchors into ortho-groups',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('anchorfiles', nargs="+", help='one of more anchor files')
    sp1.add_argument('outfile', help='output file')
    sp1.set_defaults(func = group)

    sp1 = sp.add_parser("omgprepare", 
            help = 'prepare weights file to run Sankoff OMG algorithm',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('ploidy', help='ploidy')
    sp1.add_argument('anchorfile', help='anchor file')
    sp1.add_argument('blastfile', help='blast file')
    sp1.add_argument("--norbh", action="store_true", help="Disable RBH hits")
    sp1.add_argument("--pctid", default=0, type=int, 
            help="Percent id cutoff for RBH hits")
    sp1.add_argument("--cscore", default=90, type=int,
            help="C-score cutoff for RBH hits")
    sp1.set_defaults(func = omgprepare)

    sp1 = sp.add_parser("omg", 
            help = 'generate a series of Sankoff OMG algorithm inputs',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('weightsfiles', nargs="+", help='one of more weights files')
    sp1.set_defaults(func = omg)

    sp1 = sp.add_parser("omgparse", 
            help = 'parse the OMG outputs to get gene lists',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('work', help='work')
    sp1.set_defaults(func = omgparse)

    sp1 = sp.add_parser("layout", help = 'layout the gene lists',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('omgfile', help='omgfile')
    sp1.add_argument('taxa', help='taxa')
    sp1.add_argument("--sort", help="Sort layout file based on bedfile")
    sp1.set_defaults(func = layout)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

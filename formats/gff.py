#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import os
import os.path as op
import logging
import re

from itertools import chain
from more_itertools import flatten

from jcvi.utils.cbook import AutoVivification
from jcvi.formats.base import DictFile, LineFile, must_open, is_number
from jcvi.apps.base import mkdir, parse_multi_values, need_update, sh
from jcvi.formats.gff import Gff, GffLine
from jcvi.formats.bed import Bed, BedLine
from natsort import natsorted
from jcvi.utils.range import range_minmax

valid_gene_type = 'gene'
valid_rna_type = """
    mRNA rRNA tRNA
    miRNA ncRNA lnc_RNA snRNA snoRNA pre_miRNA SRP_RNA RNase_MRP_RNA
    pseudogenic_transcript miRNA_primary_transcript
""".split()
valid_mrna_child_type = "exon CDS five_prime_UTR three_prime_UTR".split()
d1 = {x: ["match_part"] for x in 'match cDNA_match EST_match nucleotide_to_protein_match expressed_sequence_match protein_match'.split()}
d2 = {"transposable_element": ["transposon_fragment"], "transposable_element_gene": ["transposable_element"]}
d3 = {"gene": valid_rna_type, 'mRNA': valid_mrna_child_type}
d4 = {x: ['exon'] for x in valid_rna_type if x != 'mRNA'}
valid_gff_parent_child = {**d1, **d2, **d3, **d4}
valid_gff_type = set(list(valid_gff_parent_child.keys()) +
        list(chain.from_iterable(valid_gff_parent_child.values())))
valid_gff_to_gtf_type = {
    "exon": "exon", "pseudogenic_exon": "exon", "CDS": "CDS",
    "start_codon": "start_codon", "stop_codon": "stop_codon",
    "five_prime_UTR": "5UTR", "three_prime_UTR": "3UTR"
}
reserved_gff_attributes = ("ID", "Name", "Alias", "Parent", "Target",
                           "Gap", "Derives_from", "Note", "Dbxref",
                           "Ontology_term", "Is_circular")
multiple_gff_attributes = ("Parent", "Alias", "Dbxref", "Ontology_term")
safechars = " /:?~#+!$'@()*[]|"

def make_index(gff_file):
    """
    Make a sqlite database for fast retrieval of features.
    """
    import gffutils
    db_file = gff_file + ".db"

    if need_update(gff_file, db_file):
        if op.exists(db_file):
            os.remove(db_file)
        logging.debug("Indexing `{0}`".format(gff_file))
        gffutils.create_db(gff_file, db_file, force = True,
                merge_strategy="create_unique")
    else:
        logging.debug("Load index `{0}`".format(gff_file))

    return gffutils.FeatureDB(db_file)

def index(args):
    import gffutils
    gff_file = args.fi
    db_file = args.fo
    #db_file = gff_file + ".db"
    if need_update(gff_file, db_file):
        if op.exists(db_file):
            os.remove(db_file)
        logging.debug("Indexing `{0}`".format(gff_file))
        gffutils.create_db(gff_file, db_file, force = True,
                merge_strategy="create_unique")
    else:
        logging.debug("Index already made: `{0}`".format(gff_file))

def make_attributes(s, gff3=True, keep_attr_order=True):
    """
    In GFF3, the last column is typically:
    ID=cds00002;Parent=mRNA00002;
    In GFF2, the last column is typically:
    Gene 22240.t000374; Note "Carbonic anhydrase"
    """
    if gff3:
        """
        hack: temporarily replace the '+' sign in the attributes column
        with the string 'PlusSign' to prevent urlparse.parse_qsl() from
        replacing the '+' sign with a space
        """
        s = s.replace('+', 'PlusSign')
        d = parse_qs(s)
        for key in d.keys():
            d[key][0] = unquote(d[key][0].replace('PlusSign', '+').replace('"', ''))
    else:
        attributes = s.split(";")
        d = dict()
        for a in attributes:
            a = a.strip()
            if ' ' not in a:
                continue
            key, val = a.split(' ', 1)
            val = unquote(val.replace('"', '').replace('=', ' ').strip())
            if key not in d: d[key] = []
            d[key].append(val)

    for key, val in d.items():
        d[key] = list(flatten([v.split(",") for v in val]))

    return d

def addparent(args):
    """
    %prog addparent file.gff

    Merge sister features and infer parents.
    """
    p = OptionParser(addparent.__doc__)
    p.add_argument("--childfeat", default="CDS", help="Type of children feature")
    p.add_argument("--parentfeat", default="mRNA", help="Type of merged feature")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gff_file, = args
    gff = Gff(gff_file)
    data = defaultdict(list)
    for g in gff:
        if g.type != args.childfeat:
            continue
        data[g.parent].append(g)

    logging.debug("A total of {0} {1} features clustered".\
                    format(len(data), args.childfeat))

    parents = []
    for parent, dd in data.items():
        d = dd[0]
        start, end = min(x.start for x in dd), max(x.end for x in dd)
        gffline = "\t".join(str(x) for x in \
                        (d.seqid, d.source, args.parentfeat, start, end,
                         ".", d.strand, ".", "ID={0};Name={0}".format(parent)))
        parents.append(GffLine(gffline))
    parents.sort(key=lambda x: (x.seqid, x.start))
    logging.debug("Merged feature sorted")

    fw = must_open(args.outfile, "w")
    for parent in parents:
        print >> fw, parent
        parent_id = parent.id
        for d in data[parent_id]:
            if d.accn == parent_id:
                new_id = "{0}.{1}1".format(parent_id, args.childfeat)
                d.set_attr("ID", new_id)
                d.set_attr("Name", new_id, update=True)
            print >> fw, d
    fw.close()

def _fasta_slice(fasta, seqid, start, stop, strand):
    """
    Return slice of fasta, given (seqid, start, stop, strand)
    """
    _strand = 1 if strand == '+' else -1
    return fasta.sequence({'chr': seqid, 'start': start, 'stop': stop, \
        'strand': _strand})

def is_valid_codon(codon, type='start'):
    """
    Given a codon sequence, check if it is a valid start/stop codon
    """
    if len(codon) != 3:
        return False

    if type == 'start':
        if codon != 'ATG':
            return False
    elif type == 'stop':
        if not any(_codon == codon for _codon in ('TGA', 'TAG', 'TAA')):
            return False
    else:
        logging.error("`{0}` is not a valid codon type. ".format(type) + \
            "Should be one of (`start` or `stop`)")
        sys.exit()

    return True

def scan_for_valid_codon(codon_span, strand, seqid, genome, type='start'):
    """
    Given a codon span, strand and reference seqid, scan upstream/downstream
    to find a valid in-frame start/stop codon
    """
    s, e = codon_span[0], codon_span[1]
    while True:
        if (type == 'start' and strand == '+') or \
            (type == 'stop' and strand == '-'):
            s, e = s - 3, e - 3
        else:
            s, e = s + 3, e + 3

        codon = _fasta_slice(genome, seqid, s, e, strand)
        is_valid = is_valid_codon(codon, type=type)
        if not is_valid:
            if type == 'start':
                ## if we are scanning upstream for a valid start codon,
                ## stop scanning when we encounter a stop
                if is_valid_codon(codon, type='stop'):
                    return (None, None)
            elif type == 'stop':
                ## if we are scanning downstream for a valid stop codon,
                ## stop scanning when we encounter a start
                if is_valid_codon(codon, type='start'):
                    return (None, None)
            continue
        break

    return (s, e)

def cluster(args):
    """
    %prog cluster gffile

    Given a gff file of gene structures (multiple transcripts per gene locus),
    cluster/consolidate all transcripts based on shared splicing structure.

    If `slop` is enabled, clustering/consolidation will collapse any variation
    in terminal UTR lengths, keeping only the longest as representative.
    """
    from maize.utils.grouper import Grouper
    from itertools import combinations

    gffile = args.gff
    slop = args.slop
    inferUTR = args.inferUTR

    gff = make_index(gffile)

    fw = must_open(args.outfile, "w")
    print >> fw, "##gff-version	3"
    seen = {}
    for gene in gff.features_of_type('gene', order_by=('seqid', 'start')):
        g = Grouper()
        mrnas = list(combinations([mrna for mrna in gff.children(gene, featuretype='mRNA', order_by=('start'))], 2))
        if len(mrnas) > 0:
            for mrna1, mrna2 in mrnas:
                mrna1s, mrna2s = gff.children_bp(mrna1, child_featuretype='exon'), \
                    gff.children_bp(mrna2, child_featuretype='exon')
                g.join((mrna1.id, mrna1s))
                g.join((mrna2.id, mrna2s))

                if match_subfeats(mrna1, mrna2, gff, gff, featuretype='CDS'):
                    res = []
                    ftypes = ['exon'] if inferUTR else ['five_prime_UTR', 'three_prime_UTR']
                    for ftype in ftypes:
                        res.append(match_subfeats(mrna1, mrna2, gff, gff, featuretype=ftype, slop=slop))

                    if all(r == True for r in res):
                        g.join((mrna1.id, mrna1s), (mrna2.id, mrna2s))
        else:
            for mrna1 in gff.children(gene, featuretype='mRNA', order_by=('start')):
                mrna1s = gff.children_bp(mrna1, child_featuretype='exon')
                g.join((mrna1.id, mrna1s))

        print >> fw, gene
        for group in g:
            group.sort(key=lambda x: x[1], reverse=True)
            mrnas = [el[0] for el in group]
            m = mrnas[0]

            _mrnaid = []
            for x in mrnas:
                if x not in _mrnaid: _mrnaid.append(x)
            mrnaid = "{0}".format("-".join(_mrnaid))
            if mrnaid not in seen:
                seen[mrnaid] = 0
            else:
                seen[mrnaid] += 1
                mrnaid = "{0}-{1}".format(mrnaid, seen[mrnaid])

            _mrna = gff[m]
            _mrna.attributes['ID'] = [mrnaid]
            _mrna.attributes['Parent'] = [gene.id]
            children = gff.children(m, order_by='start')
            print >> fw, _mrna
            for child in children:
                child.attributes['ID'] = ["{0}".format(child.id)]
                child.attributes['Parent'] = [mrnaid]
                print >> fw, child

    fw.close()

def summary(args):
    """
    %prog summary gffile

    Print summary stats for features of different types.
    """
    from maize.formats.base import SetFile
    from maize.formats.bed import BedSummary
    from maize.utils.table import tabulate
    gff_file = args.fi
    ids = args.ids

    if ids:
        ids = SetFile(ids)
        logging.debug("Total ids loaded: {0}".format(len(ids)))

        if args.isoform:
            pids = set()
            gff = Gff(gff_file)
            for g in gff:
                if g.type != "mRNA":
                    continue
                if g.parent not in ids:
                    continue
                if "longest" not in g.attributes:
                    pids = set(x + ".1" for x in ids)
                    break
                if g.attributes["longest"][0] == "0":
                    continue
                pids.add(g.id)
            ids = pids
            logging.debug("After checking longest: {0}".format(len(ids)))

        # Collects aliases
        gff = Gff(gff_file)
        for g in gff:
            if g.name in ids:
                ids.add(g.id)
        logging.debug("Total ids including aliases: {0}".format(len(ids)))

    gff = Gff(gff_file)
    beds = defaultdict(list)
    for g in gff:
        if ids and not (g.id in ids or g.name in ids or g.parent in ids):
            continue

        beds[g.type].append(g.bedline)

    table = {}
    for type, bb in sorted(beds.items()):
        bs = BedSummary(bb)
        table[(type, "Features")] = bs.nfeats
        table[(type, "Unique bases")] = bs.unique_bases
        table[(type, "Total bases")] = bs.total_bases

    print >> sys.stdout, tabulate(table)

def orient(args):
    """
    %prog orient in.gff3 features.fasta > out.gff3

    Change the feature orientations based on translation. This script is often
    needed in fixing the strand information after mapping RNA-seq transcripts.

    You can generate the features.fasta similar to this command:

    $ %prog load --parents=EST_match --children=match_part clc.JCVIv4a.gff
    JCVI.Medtr.v4.fasta -o features.fasta
    """
    from maize.formats.fasta import longestorf
    ingff3, fastafile = args.fi, args.fasta
    idsfile = fastafile.rsplit(".", 1)[0] + ".orf.ids"
    if need_update(fastafile, idsfile):
        longestorf([fastafile, "--ids"])

    orientations = DictFile(idsfile)
    gff = Gff(ingff3)
    flipped = 0
    for g in gff:
        id = None
        for tag in ("ID", "Parent"):
            if tag in g.attributes:
                id, = g.attributes[tag]
                break
        assert id

        orientation = orientations.get(id, "+")
        if orientation == '-':
            g.strand = {"+": "-", "-": "+"}[g.strand]
            flipped += 1

        print(g)

    logging.debug("A total of {0} features flipped.".format(flipped))

def rename(args):
    """
    %prog rename in.gff3 switch.ids > reindexed.gff3

    Change the IDs within the gff3.
    """
    ingff3, switch = args.fi, args.map
    switch = DictFile(switch)

    gff = Gff(ingff3)
    for g in gff:
        id, = g.attributes["ID"]
        newname = switch.get(id, id)
        g.attributes["ID"] = [newname]

        if "Parent" in g.attributes:
            parents = g.attributes["Parent"]
            g.attributes["Parent"] = [switch.get(x, x) for x in parents]

        g.update_attributes()
        print(g)

def filter(args):
    """
    %prog filter gffile > filtered.gff

    Filter the gff file based on criteria below:
    (1) feature attribute values: [Identity, Coverage].
    You can get this type of gff by using gmap
    $ gmap -f 2 ....

    (2) Total bp length of child features
    """
    gffile = args.fi
    otype, oid, ocov = args.type, args.id, args.coverage
    cftype, clenbp = args.child_ftype, args.child_bp

    id_attr, cov_attr = "Identity", "Coverage"
    if args.nocase:
        id_attr, cov_attr = id_attr.lower(), cov_attr.lower()

    gffdb = make_index(gffile)
    bad = set()
    ptype = None
    for g in gffdb.features_of_type(otype, order_by=('seqid', 'start')):
        if not ptype:
            parent = list(gffdb.parents(g))
            ptype = parent[0].featuretype \
                if len(parent) > 0 else otype
        if cftype and clenbp:
            if gffdb.children_bp(g, child_featuretype=cftype) < clenbp:
                bad.add(g.id)
        elif oid and ocov:
            identity = float(g.attributes[id_attr][0])
            coverage = float(g.attributes[cov_attr][0])
            if identity < oid or coverage < ocov:
                bad.add(g.id)

    logging.debug("{0} bad accns marked.".format(len(bad)))
    for g in gffdb.features_of_type(ptype, order_by=('seqid', 'start')):
        if ptype != otype:
            feats = list(gffdb.children(g, featuretype=otype, order_by=('start')))
            ok_feats = [f for f in feats if f.id not in bad]
            if len(ok_feats) > 0:
                print(g)
                for feat in ok_feats:
                    print(feat)
                    for child in gffdb.children(feat, order_by=('start')):
                        print(child)
        else:
            if g.id not in bad:
                print(g)
                for child in gffdb.children(g, order_by=('start')):
                    print(child)

def fix_gsac(g, notes):
    a = g.attributes

    if g.type == "gene":
        note = a["Name"]
    elif g.type == "mRNA":
        parent = a["Parent"][0]
        note = notes[parent]
    else:
        return

    a["Name"] = a["ID"]
    a["Note"] = note
    g.update_attributes()

def fix(args):
    fi, opt = args.fi, args.opt
    gff = Gff(fi)
    if opt == 'genbank':
        for g in gff:
            if g.type == 'region':
                continue
            print(g)
    elif opt == 'tair':
        for g in gff:
            if g.type in ['protein','chromosome','transposon_fragment','transposable_element']:
                continue
            elif g.type == 'pseudogenic_transcript':
                g.type = 'mRNA'
            elif g.type == 'pseudogenic_exon':
                g.type = 'exon'
            elif g.type == 'CDS':
                parent = g.get_attr('Parent')
            print(g)
    elif opt == 'phytozome':
        gdic = dict()
        for g in gff:
            if g.type.endswith('RNA') or g.type=='gene':
                id, name = g.get_attr("ID"), g.get_attr("Name")
                #name = name.replace(".", '')
                g.set_attr("ID", name)
                gdic[id] = name
            if g.type != 'gene':
                g.set_attr("Parent", gdic[g.get_attr('Parent')])
            g.update_attributes()
            print(g)
    elif opt == 'maize':
        for g in gff:
            if g.type == 'transposable_element':
                g.type = 'transposable_element_gene'
            conf = g.get_attr('conf_class')
            if conf:
                g.set_attr("Note", "[%s]%s" % (g.get_attr("Note"), conf))
            print(g)
    elif opt == 'nam':
        for g in gff:
            if g.type in ['chromosome','scaffold']:
                continue
            if g.type == 'gene':
                g.set_attr("ID", g.get_attr("ID").replace("gene:", ""))
            elif g.type =='mRNA':
                g.set_attr("ID", g.get_attr("ID").replace("transcript:", ""))
                g.set_attr("Parent", g.get_attr("Parent").replace("gene:", ""))
            else:
                g.set_attr("Parent", g.get_attr("Parent").replace("transcript:", ""))
                if g.type == 'CDS':
                    g.set_attr("ID", '')
            g.update_attributes()
            print(g)
    elif opt == 'ensembl':
        seqtypes = ["chromosome","contig",'supercontig','biological_region','region','scaffold']
        id_map = dict()
        for g in gff:
            if g.type in seqtypes:
                continue
            elif g.type == "ncRNA_gene":
                g.type = "gene"
            elif g.type == 'pseudogene':
                g.type = 'gene'
            elif g.type == 'pseudogenic_transcript':
                valid_pseudo_biotype = { 'pseudogene':'mRNA',
                                        'tRNA_pseudogene':'tRNA' }
                biotype = g.get_attr('biotype')
                assert biotype and biotype in valid_pseudo_biotype, \
                    'unknown biotype: %s' % g.get_attr("ID")
                #g.type = valid_pseudo_biotype[biotype]
            elif g.type not in valid_gff_type:
                logging.error("type[%s] not allowed" % g.type)
                sys.exit(1)
            elif g.type in valid_mrna_child_type:
                if g.get_attr("ID"):
                    g.set_attr("ID", None)
            if g.type == 'gene':
                g.update_tag("Name", "note1")
                g.update_tag("description", "note2")
            if g.get_attr("ID"):
                ary = g.get_attr("ID").split(":")
                if len(ary) == 2:
                    g.set_attr("ID", ary[1])
            if g.get_attr("Parent"):
                ary = g.get_attr("Parent").split(":")
                if len(ary) == 2:
                    g.set_attr("Parent", ary[1])
            g.update_attributes()
            if g.type in valid_rna_type and g.get_attr("ID") == g.get_attr("Parent"):
                oid = g.get_attr("ID")
                nid = g.get_attr("ID") + ".1"
                g.set_attr("ID", nid)
                assert oid not in id_map, 'more than 2 children: %s' % oid
                id_map[oid] = nid
                g.update_attributes()
            if g.type in ['exon','CDS'] and g.get_attr("Parent") in id_map:
                g.set_attr("Parent", id_map[g.get_attr("Parent")])
                g.update_attributes()
            print(g)
    elif opt == 'mo17':
        gdic = dict()
        for g in gff:
            if g.type == 'region':
                continue
            elif g.type == 'five_P00rime_UTR':
                g.type = 'five_prime_UTR'
            elif g.type == 'three_P00rime_UTR':
                g.type = 'three_prime_UTR'
            # elif g.type == 'gene':
                # nid = g.get_attr('locus_tag').replace('Zm00014a_','Zm00014a')
                # oid = g.get_attr('ID')
                # g.set_attr('ID', nid)
                # gdic[oid] = nid
                # g.update_tag("Note", "note1")
                # g.update_tag("gene", "note2")
            # elif g.type.endswith('RNA'):
                # nid = g.get_attr('orig_transcript_id').replace('gnl|WGS:NCVQ|','').replace('Zm00014a_','Zm00014a')
                # oid = g.get_attr('ID')
                # g.set_attr('ID', nid)
                # gdic[oid] = nid
                # g.set_attr('Parent', gdic[g.get_attr('Parent')])
            # else:
                # g.set_attr('Parent', gdic[g.get_attr('Parent')])
            g.update_attributes()
            print(g)
    elif opt == 'hzs':
        for g in gff:
            if g.type == 'region':
                continue
            print(g)
    elif opt == 'w22':
        chrids = ["%d" % x for x in range(1, 11)]
        for g in gff:
            if g.type in ['chromosome', 'intron']:
                continue
            if g.seqid in chrids:
                g.seqid = "chr%s" % g.seqid
            elif g.seqid == 'chr10000001' or g.seqid == '10000001':
                g.seqid = 'unmapped'
            print(g)
    elif opt == 'ph207':
        dicc = {f"chr{x}": f"chr0{x}" for x in range(1, 10)}
        for g in gff:
            gid, par = g.get_attr("ID"), g.get_attr("Parent")
            if g.seqid in dicc:
                g.seqid = dicc[g.seqid]
            if gid:
                g.set_attr("ID", gid.replace(".v1.1", ""))
            if par:
                g.set_attr("Parent", par.replace(".v1.1", ""))
            g.update_attributes()
            print(g)
    elif opt == 'phb47':
        gdic = dict()
        for g in gff:
            if g.type.endswith('RNA'):
                gid, gname = g.get_attr("ID"), g.get_attr("Name")
                g.set_attr("ID", g.name)
                gdic[gid] = gname
            elif g.type != 'gene':
                g.set_attr("Parent", gdic[g.get_attr("Parent")])
            g.update_attributes()
            print(g)
    else:
        logging.error("unknown opt: %s" % opt)

def extract(args):
    """
    %prog extract gffile

    --contigs: Extract particular contig(s) from the gff file. If multiple contigs are
    involved, use "," to separate, e.g. "contig_12,contig_150"; or provide a file
    with multiple contig IDs, one per line
    --names: Process particular ID(s) from the gff file. If multiple IDs are
    involved, use "," to separate; or provide a file with multiple IDs, one per line
    """
    gffile = args.fi
    contigfile = args.contigs
    namesfile = args.names
    typesfile = args.types
    nametag = args.tag

    contigID = parse_multi_values(contigfile)
    names = parse_multi_values(namesfile)
    types = parse_multi_values(typesfile)
    if args.children:
        assert types is not None or names is not None, "Must set --names or --types"
        if names == None: names = list()
        populate_children(outfile, names, gffile, iter=args.children, types=types)
        return

    fp = must_open(gffile)
    for row in fp:
        atoms = row.split()
        if len(atoms) == 0:
            continue
        tag = atoms[0]
        if row[0] == "#":
            if row.strip() == "###":
                continue
            if not (tag == RegionTag and contigID and atoms[1] not in contigID):
                print >> fw, row.rstrip()
            if tag == FastaTag:
                break
            continue

        b = GffLine(row)
        attrib = b.attributes
        if contigID and tag not in contigID:
            continue
        if types and b.type in types:
            _id = b.accn
            if _id not in names:
                names.append(_id)
        if names is not None:
            if nametag not in attrib:
                continue
            if attrib[nametag][0] not in names:
                continue

        print(row.rstrip())

    if not args.fasta:
        return

    f = Fasta(gffile)
    for s in contigID:
        if s in f:
            SeqIO.write([f[s]], sys.stdout, "fasta")

def chain(args):
    """
    %prog chain gffile > chained.gff

    Fill in parent features by chaining child features and return extent of the
    child coordinates.
    """
    gffile = args.gff
    attrib_key = args.attrib_key
    attrib_list = args.attrib_list
    score_merge_op = args.score_merge_op
    break_chain = args.break_chain

    chain_ftype = args.chain_ftype
    parent_ftype = args.parent_ftype

    gffdict = {}
    fw = must_open(args.outfile, "w")
    gff = Gff(gffile)
    if break_chain:
        ctr, prev_gid = dict(), None
    for g in gff:
        if g.type != chain_ftype:
            print >> fw, g
            continue

        id = g.accn
        gid = id
        if attrib_key:
            assert attrib_key in g.attributes.keys(), \
                "Attribute `{0}` not present in GFF3".format(attrib_key)
            gid = g.get_attr(attrib_key)
        curr_gid = gid
        if break_chain:
            if prev_gid != curr_gid:
                if curr_gid not in ctr:
                    ctr[curr_gid] = 0
                else:
                    ctr[curr_gid] += 1
            gid = "{0}:{1}".format(gid, ctr[curr_gid])
        gkey = (g.seqid, gid)
        if gkey not in gffdict:
            gffdict[gkey] = { 'seqid': g.seqid,
                            'source': g.source,
                            'strand': g.strand,
                            'type': parent_ftype,
                            'coords': [],
                            'children': [],
                            'score': [],
                            'attrs': DefaultOrderedDict(set)
                          }
            gffdict[gkey]['attrs']['ID'].add(gid)

        if attrib_list:
            for a in attrib_list.split(","):
                if a in g.attributes:
                    [gffdict[gkey]['attrs'][a].add(x) for x in g.attributes[a]]
                    del g.attributes[a]

        if break_chain:
            _attrib = "Alias" if attrib_list and ("Name" not in attrib_list) else "Name"
            gffdict[gkey]['attrs'][_attrib].add(curr_gid)

        gffdict[gkey]['coords'].append((g.start, g.end))
        if score_merge_op:
            if is_number(g.score):
                gffdict[gkey]['score'].append(float(g.score))
                g.score = "."

        g.attributes["Parent"] = [gid]
        g.attributes["ID"] = ["{0}-{1}".\
                format(gid, len(gffdict[gkey]['children']) + 1)]
        g.update_attributes()
        gffdict[gkey]['children'].append(g)
        if break_chain:
            prev_gid = curr_gid

    for gkey, v in sorted(gffdict.items()):
        gseqid, key = gkey
        seqid = v['seqid']
        source = v['source']
        type = v['type']
        strand = v['strand']
        start, stop = range_minmax(gffdict[gkey]['coords'])

        score = "."
        if score_merge_op:
            v['score'].sort()
            if score_merge_op == "sum":
                score = sum(v['score'])
            elif score_merge_op == "min":
                score = min(v['score'])
            elif score_merge_op == "max":
                score = max(v['score'])
            elif score_merge_op == "mean":
                score = sum(v['score'], 0.0)/len(v['score'])
            elif score_merge_op == "collapse":
                score = ",".join((str(x) for x in v['score']))

        g = GffLine("\t".join(str(x) for x in [seqid, source, type, start, stop, \
            score, strand, ".", None]))
        g.attributes = v['attrs']
        g.update_attributes()

        print >> fw, g

        for child in gffdict[gkey]['children']:
            print >> fw, child

    fw.close()

def format(args):
    """
    %prog format gffile > formatted.gff

    Read in the gff and print it out, changing seqid, etc.
    """
    from maize.formats.obo import load_GODag, validate_term

    gffile = args.gff
    mapfile = args.seqid
    names = args.name
    note = args.note
    source = args.source
    ftype = args.type
    attrib_files = args.attrib_files.split(",") if args.attrib_files else None
    dbxref_files = args.dbxref_files.split(",") if args.dbxref_files else None
    remove_attrs = args.remove_attrs.split(",") if args.remove_attrs else None
    process_ftype = args.process_ftype.split(",") if args.process_ftype else None
    gsac = args.gsac
    assert not (args.unique and args.duptype), \
        "Cannot use `--unique` and `--chaindup` together"
    assert not(args.type and args.duptype), \
        "Cannot use `--type` and `--chaindup` together"
    unique = args.unique
    duptype = args.duptype
    fixphase = args.fixphase
    phaseT = {"1":"2", "2":"1"}
    remove_feats = args.remove_feats.split(",") if args.remove_feats else None
    remove_feats_by_ID = None
    if args.remove_feats_by_ID:
        remove_feats_by_ID = LineFile(args.remove_feats_by_ID, load=True).lines \
                if op.isfile(args.remove_feats_by_ID) else \
                args.remove_feats_by_ID.split(",")
    strict = False if args.nostrict else True
    make_gff_store = True if gffile in ("-", "stdin") else args.make_gff_store
    assert not (args.copy_id_attr_to_name and args.invent_name_attr), \
        "Cannot use `--copy_id_attr_to_name` and `--invent_name_attr` together"
    copy_id_attr_to_name = args.copy_id_attr_to_name
    invent_name_attr = args.invent_name_attr
    invent_protein_feat = args.invent_protein_feat
    compute_signature = False

    outfile = args.outfile

    mapping = None
    mod_attrs = set()
    if mapfile and op.isfile(mapfile):
        mapping = DictFile(mapfile, delimiter="\t", strict=strict)
        mod_attrs.add("ID")
    if note:
        note = DictFile(note, delimiter="\t", strict=strict)
        mod_attrs.add("Note")
    if source and op.isfile(source):
        source = DictFile(source, delimiter="\t", strict=strict)
    if ftype and op.isfile(ftype):
        ftype = DictFile(ftype, delimiter="\t", strict=strict)
    if names:
        names = DictFile(names, delimiter="\t", strict=strict)
        mod_attrs.add("Name")

    if attrib_files:
        attr_values = {}
        for fn in attrib_files:
            attr_name = op.basename(fn).rsplit(".", 1)[0]
            if attr_name not in reserved_gff_attributes:
                attr_name = attr_name.lower()
            attr_values[attr_name] = DictFile(fn, delimiter="\t", strict=strict)
            mod_attrs.add(attr_name)
    if dbxref_files:
        dbxref_values = {}
        for fn in dbxref_files:
            dbtag = op.basename(fn).rsplit(".", 1)[0]
            dbxref_values[dbtag] = DictFile(fn, delimiter="\t", strict=strict)
        mod_attrs.add("Dbxref")

    if remove_attrs:
        mod_remove_attrs = []
        for remove_attr in remove_attrs:
            if remove_attr in mod_attrs:
                mod_remove_attrs.append(remove_attr)

        if mod_remove_attrs:
            logging.error("Attributes `{0}` cannot be removed and modified".format( \
                    ",".join(mod_remove_attrs)))
            sys.exit()

    if gsac:  # setting gsac will force IDs to be unique
        unique = True
        notes = {}

    remove = set()
    if unique or duptype or remove_feats or remove_feats_by_ID \
            or args.multiparents == "merge" or invent_name_attr or make_gff_store \
            or invent_protein_feat:
        if unique:
            dupcounts = defaultdict(int)
            seen = defaultdict(int)
            newparentid = {}
        elif duptype:
            dupranges = AutoVivification()
            skip = defaultdict(int)
        if args.multiparents == "merge":
            merge_feats = AutoVivification()
        if invent_name_attr:
            ft = GffFeatureTracker()
        elif copy_id_attr_to_name:
            pass
        if invent_protein_feat:
            cds_track = {}
        if args.multiparents == "merge" or invent_name_attr:
            make_gff_store = compute_signature = True
        gff = Gff(gffile, keep_attr_order=(not args.no_keep_attr_order), \
                make_gff_store=make_gff_store, compute_signature=compute_signature, \
                strict=strict)
        for g in gff:
            if process_ftype and g.type not in process_ftype:
                continue
            id = g.accn
            if remove_feats and g.type in remove_feats:
                remove.add(id)
            if remove_feats_by_ID and id in remove_feats_by_ID:
                remove.add(id)
            if unique:
                dupcounts[id] += 1
            elif duptype and g.type == duptype:
                dupranges[g.seqid][id][g.idx] = (g.start, g.end)
            if args.multiparents == "merge" and g.type != "CDS": #don't merge CDS
                pp = g.get_attr("Parent", first=False)
                if pp and len(pp) > 0:
                    for parent in pp:
                        if parent not in remove:
                            sig = g.sign
                            if sig not in merge_feats:
                                merge_feats[sig]['parents'] = []
                            if parent not in merge_feats[sig]['parents']:
                                merge_feats[sig]['parents'].append(parent)
            if invent_name_attr:
                parent, iso = atg_name(g.get_attr("Parent"), retval="locus,iso")
                if not parent:
                    parent = g.get_attr("Parent")
                ft.track(parent, g)
            if invent_protein_feat:
                if g.type == 'CDS':
                    cds_parent = g.get_attr("Parent")
                    if cds_parent not in cds_track:
                        cds_track[cds_parent] = []
                    cds_track[cds_parent].append((g.start, g.end))

    if args.verifySO:
        so = load_GODag()
        valid_soterm = {}

    fw = must_open(outfile, "w")
    if not make_gff_store:
        gff = Gff(gffile, keep_attr_order=(not args.no_keep_attr_order), \
                strict=strict)
    for g in gff:
        if process_ftype and g.type not in process_ftype:
            print >> fw, g
            continue

        id = g.accn

        if args.multiparents == "merge" and g.type != "CDS": #don't merge CDS
            sig = g.sign
            if len(merge_feats[sig]['parents']) > 1:
                if 'candidate' not in merge_feats[sig]:
                    merge_feats[sig]['candidate'] = id
                    g.set_attr("Parent", merge_feats[sig]['parents'])
                else:
                    continue

        if remove_feats or remove_feats_by_ID:
            if id in remove:
                continue
            else:
                if "Parent" in g.attributes:
                    keep, parent = [], g.get_attr("Parent", first=False)
                    for i, pid in enumerate(parent):
                        if pid not in remove:
                            keep.append(parent[i])
                        else:
                            remove.add(id)
                    if len(keep) == 0:
                        continue
                    parent = g.set_attr("Parent", keep)

        if remove_attrs:
            for remove_attr in remove_attrs:
                if remove_attr in g.attributes:
                    g.set_attr(remove_attr, None)

        if args.verifySO:
            if g.type not in valid_soterm:
                valid_soterm[g.type] = validate_term(g.type, so=so, method=args.verifySO)
            ntype = valid_soterm[g.type]
            if ntype and g.type != ntype:
                g.type = ntype

        origid = g.seqid
        if fixphase:
            phase = g.phase
            g.phase = phaseT.get(phase, phase)

        if mapfile:
            if isinstance(mapping, dict):
                if origid in mapping:
                    g.seqid = mapping[origid]
                else:
                    logging.error("{0} not found in `{1}`. ID unchanged.".\
                            format(origid, mapfile))
            else:
                g.seqid = mapfile

        if source:
            if isinstance(source, dict) and g.source in source:
                g.source = source[g.source]
            else:
                g.source = source

        if names:
            if id in names:
                g.set_attr("Name", names[id])

        if note:
            name = g.get_attr("Name")
            tag = None
            if id in note:
                tag = note[id]
            elif name and name in note:
                tag = note[name]

            if tag:
                g.set_attr("Note", tag, update=False)

        if attrib_files:
            for attr_name in attr_values:
                name = g.get_attr("Name")
                if id in attr_values[attr_name]:
                    g.set_attr(attr_name, attr_values[attr_name][id])
                elif name and name in attr_values[attr_name]:
                    g.set_attr(attr_name, attr_values[attr_name][name])

        if dbxref_files:
            for dbtag in dbxref_values:
                if id in dbxref_values[dbtag]:
                    g.set_attr("Dbxref", dbxref_values[dbtag][id], dbtag=dbtag, append=True)

        if unique:
            if dupcounts[id] > 1:
                seen[id] += 1
                old_id = id
                id = "{0}-{1}".format(old_id, seen[old_id])
                newparentid[old_id] = id
                g.set_attr("ID", id)

            if "Parent" in g.attributes:
                parent = g.attributes["Parent"][0]
                if dupcounts[parent] > 1:
                    g.set_attr("Parent", newparentid[parent])

        if duptype:
            if duptype == g.type and len(dupranges[g.seqid][id]) > 1:
                p = sorted(dupranges[g.seqid][id])
                s, e = dupranges[g.seqid][id][p[0]][0:2]  # get coords of first encountered feature
                if g.start == s and g.end == e and p[0] == g.idx:
                    r = [dupranges[g.seqid][id][x] for x in dupranges[g.seqid][id]]
                    g.start, g.end = range_minmax(r)
                else:
                    skip[(g.seqid, g.idx, id, g.start, g.end)] = 1

        if gsac and g.type == "gene":
            notes[id] = g.attributes["Name"]

        if ftype:
            if isinstance(ftype, dict) and g.type in ftype:
                g.type = ftype[g.type]
            else:
                g.type = ftype

        if invent_name_attr:
            ft.store_symbol(g)
            if re.search(ft.ftype, g.type):
                parent, iso = atg_name(g.get_attr("Parent"), retval="locus,iso")
                if not parent:
                    parent = g.get_attr("Parent")
                if parent in ft.tracker:
                    fidx = ft.feat_index(parent, g.type, g.strand, (g.start, g.end, g.sign))
                    symbol = ft.get_symbol(parent)
                    attr = "ID" if symbol == parent else "Name"
                    g.set_attr(attr, "{0}:{1}:{2}".format(symbol, g.type, fidx + 1))
                    if args.multiparents == "merge" and attr == "Name":
                        g.set_attr("ID", "{0}:{1}:{2}".format(parent, g.type, fidx + 1))
        elif copy_id_attr_to_name:
            if "Name" not in g.attributes.keys():
                g.set_attr("Name", g.get_attr("ID"))

        protein_feat = None
        if invent_protein_feat:
            if g.type == 'mRNA':
                if id in cds_track:
                    pstart, pstop = range_minmax(cds_track[id])
                    protein_feat = GffLine("\t".join(str(x) for x in [g.seqid, g.source, "protein", pstart, pstop, \
                            ".", g.strand, ".", "ID={0}-Protein;Name={0};Derives_from={0}".format(id)]))
            elif g.type == 'CDS':
                parent = g.get_attr("Parent")
                if parent in cds_track:
                    _parent = [parent, "{0}-Protein".format(parent)]
                    g.set_attr("Parent", _parent)

        pp = g.get_attr("Parent", first=False)
        if args.multiparents == "split" and (pp and len(pp) > 1) and g.type != "CDS":  # separate features with multiple parents
            id = g.get_attr("ID")
            for i, parent in enumerate(pp):
                if id: g.set_attr("ID", "{0}-{1}".format(id, i + 1))
                g.set_attr("Parent", parent, update=True, urlquote=True)
                if gsac:
                    fix_gsac(g, notes)
                print >> fw, g
        else:
            if g.gff3 and not args.gff3:
                args.gff3 = True
            g.update_attributes(gff3=args.gff3)
            if gsac:
                fix_gsac(g, notes)
            if duptype == g.type and skip[(g.seqid, g.idx, id, g.start, g.end)] == 1:
                continue
            print >> fw, g
            if g.type == 'mRNA' and invent_protein_feat:
                print >> fw, protein_feat

    fw.close()

def fixboundaries(args):
    """
    %prog fixboundaries gffile --type="gene" --child_ftype="mRNA" > gffile.fixed

    Adjust the boundary coordinates of parents features based on
    range chained child features, extracting their min and max values
    """
    gffile = args.fi
    gffdb = make_index(gffile)

    for f in gffdb.all_features(order_by=('seqid', 'start')):
        if f.featuretype == args.type:
            child_coords = []
            for cftype in args.child_ftype.split(","):
                for c in gffdb.children(f, featuretype=cftype, order_by=('start')):
                    child_coords.append((c.start, c.stop))
            f.start, f.stop = range_minmax(child_coords)
        print(f)

def get_piles(allgenes):
    """
    Before running uniq, we need to compute all the piles. The piles are a set
    of redundant features we want to get rid of. Input are a list of GffLines
    features. Output are list of list of features distinct "piles".
    """
    from maize.utils.range import Range, range_piles

    ranges = [Range(a.seqid, a.start, a.end, 0, i) \
                    for i, a in enumerate(allgenes)]

    for pile in range_piles(ranges):
        yield [allgenes[x] for x in pile]

def match_span(f1, f2):
    return (f1.start == f2.start) and (f1.stop == f2.stop)

def match_ftype(f1, f2):
    return f1.featuretype == f2.featuretype

def match_nchildren(f1c, f2c):
    return len(f1c) == len(f2c)

def match_child_ftype(f1c, f2c):
    from collections import Counter

    return len(set(Counter(i.featuretype for i in f1c).keys()) ^ \
            set(Counter(i.featuretype for i in f2c).keys()))

def match_Nth_child(f1c, f2c, N=1, slop=False):
    i = N - 1
    f1, f2 = f1c[i], f2c[i]

    if slop:
        if 1 == len(f1c):
            if f1.featuretype.endswith('UTR'):
                if f1.strand == '+':
                    Npos = "F" if f1.featuretype.startswith('five_prime') else "L"
                elif f1.strand == '-':
                    Npos = "L" if f1.featuretype.startswith('five_prime') else "F"
            elif f1.featuretype == 'exon':
                return not match_span(f1, f2)
        elif N == 1: Npos = "F"
        elif N == len(f1c): Npos = "L"

        if Npos == "F":
            return f1.stop == f2.stop
        elif Npos == "L":
            return f1.start == f2.start

    return match_span(f1, f2)

def match_subfeats(f1, f2, dbx1, dbx2, featuretype=None, slop=False):
    """
    Given 2 gffutils features located in 2 separate gffutils databases,
    iterate through all subfeatures of a certain type and check whether
    they are identical or not

    The `slop` parameter allows for variation in the terminal UTR region
    """
    f1c, f2c = list(dbx1.children(f1, featuretype=featuretype, order_by='start')), \
            list(dbx2.children(f2, featuretype=featuretype, order_by='start'))

    lf1c, lf2c = len(f1c), len(f2c)
    if match_nchildren(f1c, f2c):
        if lf1c > 0 and lf2c > 0:
            exclN = set()
            if featuretype.endswith('UTR') or featuretype == 'exon':
                N = []
                if featuretype.startswith('five_prime'):
                    N = [1] if f1.strand == "+" else [lf1c]
                elif featuretype.startswith('three_prime'):
                    N = [lf1c] if f1.strand == "+" else [1]
                else:   # infer UTR from exon collection
                    N = [1] if 1 == lf1c else [1, lf1c]

                for n in N:
                    if match_Nth_child(f1c, f2c, N=n, slop=slop):
                        exclN.add(n-1)
                    else:
                        return False

            for i, (cf1, cf2) in enumerate(zip(f1c, f2c)):
                if i in exclN: continue
                if not match_span(cf1, cf2):
                    return False
    else:
        if (lf1c, lf2c) in [(0, 1), (1, 0)] and slop \
                and featuretype.endswith('UTR'):
            return True

        return False

    return True

def import_feats(gffile, type="gene"):
    gff = Gff(gffile)
    allgenes = []
    for g in gff:
        if g.type != type:
            continue
        allgenes.append(g)

    logging.debug("A total of {0} {1} features imported.".format(len(allgenes), type))
    allgenes.sort(key=lambda x: (x.seqid, x.start))
    return allgenes

def populate_children(outfile, ids, gffile, iter="2", types=None):
    ids = set(ids)
    fw = must_open(outfile, "w")
    logging.debug("A total of {0} features selected.".format(len(ids)))
    logging.debug("Populate children. Iteration 1..")
    gff = Gff(gffile)
    children = set()
    for g in gff:
        if types and g.type in types:
            _id = g.accn
            if _id not in ids:
                ids.add(_id)
        if "Parent" not in g.attributes:
            continue
        for parent in g.attributes["Parent"]:
            if parent in ids:
                children.add(g.accn)

    if iter == "2":
        logging.debug("Populate grand children. Iteration 2..")
        gff = Gff(gffile)
        for g in gff:
            if "Parent" not in g.attributes:
                continue
            for parent in g.attributes["Parent"]:
                if parent in children:
                    children.add(g.accn)

    logging.debug("Populate parents..")
    gff = Gff(gffile)
    parents = set()
    for g in gff:
        if g.accn not in ids:
            continue
        if "Parent" not in g.attributes:
            continue
        for parent in g.attributes["Parent"]:
            parents.add(parent)

    combined = ids | children | parents
    logging.debug("Original: {0}".format(len(ids)))
    logging.debug("Children: {0}".format(len(children)))
    logging.debug("Parents: {0}".format(len(parents)))
    logging.debug("Combined: {0}".format(len(combined)))

    logging.debug("Filter gff file..")
    gff = Gff(gffile)
    seen = set()
    for g in gff:
        accn = g.accn
        if accn in seen:
            continue
        if accn in combined:
            seen.add(accn)
            print >> fw, g
    fw.close()

def get_upstream_coords(uSite, uLen, seqlen, feat, children_list, gffdb):
    """
    Subroutine takes upstream site, length, reference sequence length,
    parent mRNA feature (GffLine object), list of child feature types
    and a GFFutils.GFFDB object as the input

    If upstream of TSS is requested, use the parent feature coords
    to extract the upstream sequence

    If upstream of TrSS is requested,  iterates through all the
    children (CDS features stored in the sqlite GFFDB) and use child
    feature coords to extract the upstream sequence

    If success, returns the upstream start and stop coordinates
    else, returns None
    """
    if uSite == "TSS":
        (upstream_start, upstream_stop) = \
                (feat.start - uLen, feat.start - 1) \
                if feat.strand == "+" else \
                (feat.end + 1, feat.end + uLen)
    elif uSite == "TrSS":
        children = []
        for c in gffdb.children(feat.id, 1):

            if c.featuretype not in children_list:
                continue
            children.append((c.start, c.stop))

        if not children:
            logging.error("%s has no children with type %s" % (feat.id, ','.join(children_list)))
            return None, None

        cds_start, cds_stop = range_minmax(children)
        (upstream_start, upstream_stop) = \
                (cds_start - uLen, cds_start - 1) \
                if feat.strand == "+" else \
                (cds_stop + 1, cds_stop + uLen)

    if feat.strand == "+" and upstream_start < 1:
        upstream_start = 1
    elif feat.strand == "-" and upstream_stop > seqlen:
        upstream_stop = seqlen

    actual_uLen = upstream_stop - upstream_start + 1
    if actual_uLen < uLen:
        logging.error("sequence upstream of {0} ({1} bp) is less than upstream length {2}" .format(feat.id, actual_uLen, uLen))
        return None, None

    return upstream_start, upstream_stop

def parse_feature_param(feature):
    """
    Take the --feature param (coming from gff.load() and parse it.
    Returns feature, parents and children terms.

    Also returns length of upstream sequence (and start site) requested

    If erroneous, returns a flag and error message to be displayed on exit
    """
    # can request upstream sequence only from the following valid sites
    valid_upstream_sites = ["TSS", "TrSS"]

    upstream_site, upstream_len = None, None
    flag, error_msg = None, None
    parents, children = None, None
    if re.match(r'upstream', feature):
        parents, children = "mRNA", "CDS"
        feature, upstream_site, upstream_len = re.search(r'([A-z]+):([A-z]+):(\S+)', \
                feature).groups()

        if not is_number(upstream_len):
            flag, error_msg = 1, "Error: upstream len `" + upstream_len + "` should be an integer"

        upstream_len = int(upstream_len)
        if(upstream_len < 0):
            flag, error_msg = 1, "Error: upstream len `" + str(upstream_len) + "` should be > 0"

        if not upstream_site in valid_upstream_sites:
            flag, error_msg = 1, "Error: upstream site `" + upstream_site + "` not valid." + \
                    " Please choose from " + valid_upstream_sites
    elif feature == "CDS":
        parents, children = "mRNA", "CDS"
    else:
        flag, error_msg = 1, "Error: unrecognized argument --feature=" + feature

    return feature, parents, children, upstream_site, upstream_len, flag, error_msg

def pick_longest(args):
    gffile = args.fi
    g = make_index(gffile)
    for gene in g.features_of_type("gene"):
        rnaP, rnalenP = None, 0
        for rna in g.children(gene, level = 1):
            #if gene['ID'][0].startswith("EPlORYSAT000373610"):
            #    logging.debug(gene['ID'][0], rna.featuretype, rna['ID'])
            rnalen = 0
            if rna.featuretype == "mRNA":
                for cds in g.children(rna, featuretype = 'CDS'):
                    rnalen += cds.end - cds.start + 1
            else:
                for exon in g.children(rna, featuretype = 'exon'):
                    rnalen += exon.end - exon.start + 1
            if rnalen > rnalenP:
                rnaP, rnalenP = rna, rnalen
        if rnalenP == 0:
            logging.debug("%s skipped" % gene.attributes["ID"])
            continue
        print(gene)
        print(rnaP)
        for child in g.children(rnaP, order_by=('start')):
            print(child)

def gff2gtf(args):
    """
    %prog gtf gffile

    Convert gff to gtf file. In gtf, only exon/CDS features are important. The
    first 8 columns are the same as gff, but in the attributes field, we need to
    specify "gene_id" and "transcript_id".
    """
    gff = Gff(args.fi)
    transcript_info = AutoVivification()
    for g in gff:
        if g.type.endswith(("RNA", "transcript")):
            if "ID" in g.attributes and "Parent" in g.attributes:
                transcript_id = g.get_attr("ID")
                gene_id = g.get_attr("Parent")
            elif "mRNA" in g.attributes and "Gene" in g.attributes:
                transcript_id = g.get_attr("mRNA")
                gene_id = g.get_attr("Gene")
            else:
                transcript_id = g.get_attr("ID")
                gene_id = transcript_id
            transcript_info[transcript_id]["gene_id"] = gene_id
            transcript_info[transcript_id]["gene_type"] = g.type
            continue

        if g.type not in valid_gff_to_gtf_type.keys():
            continue

        try:
            transcript_id = g.get_attr("Parent", first=False)
        except IndexError:
            transcript_id = g.get_attr("mRNA", first=False)

        g.type = valid_gff_to_gtf_type[g.type]
        for tid in transcript_id:
            if tid not in transcript_info: continue
            gene_type = transcript_info[tid]["gene_type"]
            if not gene_type.endswith("RNA") and not gene_type.endswith("transcript"):
                continue
            gene_id = transcript_info[tid]["gene_id"]
            biotype = transcript_info[tid]["gene_type"]
            g.attributes = dict(gene_id=[gene_id], transcript_id=[tid], gene_biotype=[biotype])
            g.update_attributes(gtf=True, urlquote=False)

            print(g)

def gff2tsv(args):
    gffile = args.fi
    db = make_index(gffile)
    print("\t".join('gid tid ttype etype chrom start end srd'.split(" ")))
    for gene in db.features_of_type('gene'):
        gid = gene.id
        chrom, srd = gene.chrom, gene.strand
        for rna in db.children(gene, 1):
            tid, ttype = rna.id, rna.featuretype
            for feat in db.children(rna):
                chrom0, start, end, srd0, etype = feat.chrom, feat.start, \
                        feat.end, feat.strand, feat.featuretype
                assert chrom0 == chrom, 'sub-feat chrom error: %s' % rna["ID"]
                if srd0 != srd:
                    logging.warning("sub-feat strand != parent: %s" % rna["ID"])
                ary = (gid, tid, ttype, etype, chrom, start, end, srd0)
                print("\t".join(str(x) for x in ary))

def gff2gb(args):
    """
    %prog gb gffile fastafile

    Convert GFF3 to Genbank format. Recipe taken from:
    <http://www.biostars.org/p/2492/>
    """
    from Bio.Alphabet import generic_dna
    try:
        from BCBio import GFF
    except ImportError:
        print >> sys.stderr, "You need to install dep first: $ easy_install bcbio-gff"
    gff_file, fasta_file = args.fi, args.fasta
    pf = op.splitext(gff_file)[0]
    out_file = pf + ".gb"
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(gff_iter, out_file, "genbank")

def gff2bed12(args):
    """
    %prog bed12 gffile > bedfile

    Produce bed12 file for coding features. The exons will be converted to blocks.
    The CDS range will be shown between thickStart to thickEnd. For reference,
    bed format consists of the following fields:

    1. chrom
    2. chromStart
    3. chromEnd
    4. name
    5. score
    6. strand
    7. thickStart
    8. thickEnd
    9. itemRgb
    10. blockCount
    11. blockSizes
    12. blockStarts
    """
    gffile = args.fi
    parent, block, thick = args.parent, args.block, args.thick

    g = make_index(gffile)
    for f in g.features_of_type(parent):
        chrom = f.chrom
        chromStart = f.start - 1
        chromEnd = f.stop
        name = f.id
        score = 0
        strand = f.strand
        thickStart = chromStart
        thickEnd = chromEnd
        blocks = []

        for c in g.children(name, 1):
            cstart, cend = c.start - 1, c.stop
            if c.featuretype == block:
                blockStart = cstart - chromStart
                blockSize = cend - cstart
                blocks.append((blockStart, blockSize))
            elif c.featuretype == thick:
                thickStart = min(thickStart, cstart)
                thickEnd = max(thickEnd, cend)

        if len(blocks) == 0:
            continue
        blocks.sort()
        blockStarts, blockSizes = zip(*blocks)
        blockCount = len(blocks)
        blockSizes = ",".join(str(x) for x in blockSizes) + ","
        blockStarts = ",".join(str(x) for x in blockStarts) + ","
        itemRgb = 0

        print("\t".join(str(x) for x in (chrom, chromStart, chromEnd, \
                name, score, strand, thickStart, thickEnd, itemRgb,
                blockCount, blockSizes, blockStarts)))

def gff2fas(args):
    '''
    %prog load gff_file fasta_file [--arguments]

    Parses the selected features out of GFF, with subfeatures concatenated.
    For example, to get the CDS sequences, do this:
    $ %prog load athaliana.gff athaliana.fa --parents mRNA --children CDS

    To get 500bp upstream of a genes Transcription Start Site (TSS), do this:
    $ %prog load athaliana.gff athaliana.fa --feature=upstream:TSS:500

    Switch TSS with TrSS for Translation Start Site.
    '''
    from datetime import datetime as dt
    from pyfaidx import Fasta
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    # can request output fasta sequence id to be picked from following attributes
    valid_id_attributes = ["ID", "Name", "Parent", "Alias", "Target"]
    gff_file, fasta_file = args.gff, args.fasta

    if args.feature:
        args.feature, args.parent, args.children, upstream_site, upstream_len, \
                flag, error_msg = parse_feature_param(args.feature)
        if flag:
            sys.exit(error_msg)

    parents = set(args.parents.split(','))
    children_list = set(args.children.split(','))

    """
    In a situation where we want to extract sequence for only the top-level
    parent feature, specify feature type of parent == child
    """
    skipChildren = True if len(parents.symmetric_difference(children_list)) == 0 \
            else False

    id_attr = args.id_attribute
    desc_attr = args.desc_attribute
    sep = args.sep

    import gffutils
    g = make_index(gff_file)
    f = Fasta(fasta_file)
    seqlen = {}
    for seqid in f.keys():
        seqlen[seqid] = len(f[seqid])

    for feat in g.features_of_type(parents):
        desc = ""
        if desc_attr:
            fparent = feat.attributes['Parent'][0] \
                if 'Parent' in feat.attributes else None
            if fparent:
                try:
                    g_fparent = g[fparent]
                except gffutils.exceptions.FeatureNotFoundError:
                    logging.error("{} not found in index .. skipped".format(fparent))
                    continue
                if desc_attr in g_fparent.attributes:
                    desc = ",".join(g_fparent.attributes[desc_attr])
            elif desc_attr in feat.attributes:
                desc = ",".join(feat.attributes[desc_attr])

        if args.full_header:
            desc_parts = []
            desc_parts.append(desc)

            if args.conf_class and 'conf_class' in feat.attributes:
                desc_parts.append(feat.attributes['conf_class'][0])

            if args.full_header == "tair":
                orient = "REVERSE" if feat.strand == "-" else "FORWARD"
                feat_coords = "{0}:{1}-{2} {3} LENGTH=[LEN]".format(feat.seqid, \
                    feat.start, feat.end, orient)
            else:
                (s, e) = (feat.start, feat.end) if (feat.strand == "+") \
                        else (feat.end, feat.start)
                feat_coords = "{0}:{1}-{2}".format(feat.seqid, s, e)
            desc_parts.append(feat_coords)

            datestamp = args.datestamp if args.datestamp else \
                    "{0}{1}{2}".format(dt.now().year, dt.now().month, dt.now().day)
            desc_parts.append(datestamp)

            desc = sep.join(str(x) for x in desc_parts)
            desc = "".join(str(x) for x in (sep, desc)).strip()

        if args.feature == "upstream":
            upstream_start, upstream_stop = get_upstream_coords(upstream_site, upstream_len, \
                     seqlen[feat.seqid], feat, children_list, g)

            if not upstream_start or not upstream_stop:
                continue

            rc = True if feat.strand == '-' else False
            feat_seq = f.get_seq(feat.seqid, upstream_start, upstream_stop, rc = rc).seq

            (s, e) = (upstream_start, upstream_stop) \
                    if feat.strand == "+" else \
                     (upstream_stop, upstream_start)
            upstream_seq_loc = str(feat.seqid) + ":" + str(s) + "-" + str(e)
            desc = sep.join(str(x) for x in (desc, upstream_seq_loc, \
                    "FLANKLEN=" + str(upstream_len)))
        else:
            children = []
            if not skipChildren:
                for c in g.children(feat.id, 1):
                    if c.featuretype not in children_list:
                        continue
                    rc = True if c.strand == '-' else False
                    child = f.get_seq(c.seqid, c.start, c.stop, rc = rc).seq
                    children.append((child, c))

                if not children:
                    logging.error("%s has no children with type %s" % (feat.id, ','.join(children_list)))
                    continue
            else:
                rc = True if feat.strand == '-' else False
                child = f.get_seq(feat.seqid, feat.start, feat.stop, rc = rc).seq
                children.append((child, feat))

            # sort children in incremental position
            children.sort(key=lambda x: x[1].start)
            # reverse children if negative strand
            if feat.strand == '-':
                children.reverse()
            feat_seq = ''.join(x[0] for x in children)

        desc = desc.replace("\"", "")

        id = ",".join(feat.attributes[id_attr]) if id_attr \
                and feat.attributes[id_attr] else \
                feat.id

        if args.full_header == "tair":
            desc = desc.replace("[LEN]", str(len(feat_seq)))

        rcd = SeqRecord(Seq(feat_seq), id=id, description=desc)
        SeqIO.write([rcd], sys.stdout, 'fasta')

def frombed(args):
    """
    %prog frombed bed_file [--arguments] > gff_file

    Convert bed to gff file. In bed, the accn will convert to key='ID'
    Default type will be `match` and default source will be `source`
    """
    p = OptionParser(frombed.__doc__)
    p.add_argument("--type", default="match",
                 help="GFF feature type [default: %default]")
    p.add_argument("--source", default="default",
                help="GFF source qualifier [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    bed = Bed(bedfile)

    for b in bed:
        print(b.gffline(type=args.type, source=args.source))

def fromgtf(args):
    """
    %prog fromgtf gtffile

    Convert gtf to gff file. In gtf, the "transcript_id" will convert to "ID=",
    the "transcript_id" in exon/CDS feature will be converted to "Parent=".
    """
    gtffile = args.fi
    gff = Gff(gtffile)
    transcript_id, gene_id = args.transcript_id, args.gene_id
    nfeats = 0
    for g in gff:
        if g.type in ("transcript", "mRNA"):
            g.type = "mRNA"
            g.update_tag(transcript_id, "ID")
            g.update_tag("mRNA", "ID")
            g.update_tag(gene_id, "Parent")
            g.update_tag("Gene", "Parent")
        elif g.type in ("exon", "CDS") or "UTR" in g.type:
            g.update_tag("transcript_id", "Parent")
            g.update_tag(g.type, "Parent")
        elif g.type == "gene":
            g.update_tag(gene_id, "ID")
            g.update_tag("Gene", "ID")
        elif g.type in ['start_codon', 'stop_codon']:
            continue
        else:
            assert 0, "Don't know how to deal with {0}".format(g.type)

        g.update_attributes()
        print(g)
        nfeats += 1

    logging.debug("A total of {0} features written.".format(nfeats))

def merge(args):
    """
    %prog merge gffiles

    Merge several gff files into one. When only one file is given, it is assumed
    to be a file with a list of gff files.
    """
    print("##gff-version 3")
    attrs = ['ID', 'Parent', 'gene_id']
    for line in must_open(args.config):
        line = line.strip(" \t\n\r")
        if line == "":
            continue
        (pre, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            logging.error("no access to input file: %s" % fi)
            sys.exit(1)
        fhi = must_open(fi)
        for ln in fhi:
            if ln.startswith("#"):
                continue
            gl = GffLine(ln)
            gl.seqid = "%s|%s" % (pre, gl.seqid)
            for attr in attrs:
                v = gl.get_attr(attr)
                if v != None:
                    gl.set_attr(attr, "%s|%s" % (pre, v), update=True)
            print(gl)

def note(args):
    """
    %prog note gffile > tabfile

    Extract certain attribute field for each feature.
    """
    gffile = args.fi
    type = args.type
    if type:
        type = type.split(",")

    g = make_index(gffile)
    exoncounts = {}
    if args.exoncount:
        for feat in g.features_of_type("mRNA"):
            nexons = 0
            for c in g.children(feat.id, 1):
                if c.featuretype != "exon":
                    continue
                nexons += 1
            exoncounts[feat.id] = nexons

    attrib = args.attribute.split(",")

    gff = Gff(gffile)
    seen = set()
    AED = args.AED
    print("\t".join(['id'] + attrib))
    for g in gff:
        if type and g.type not in type:
            continue
        if AED is not None and float(g.attributes["_AED"][0]) > AED:
            continue
        keyval = [g.accn]
        for x in attrib:
            if x in g.attributes:
                keyval.append(",".join(g.attributes[x]))
            else:
                keyval.append('')
        if exoncounts:
            nexons = exoncounts.get(g.accn, 0)
            keyval.append(str(nexons))
        keyval = tuple(keyval)
        if keyval not in seen:
            print("\t".join(keyval))
            seen.add(keyval)

def splicecov(args):
    """
    %prog splicecov annotation.gff3 junctions.bed

    Given an annotation GFF file (containing introns) and a
    TopHat junctions.bed file (preprocessed using formats.bed.juncs(),
    each intron gets tagged with the JUNC identifier and read coverage.

    Output is a summary table listing for each gene locus, the isoform number,
    number of splice junctions and {average, median, min & max} read coverage
    across the junctions.
    """
    from tempfile import mkstemp
    from pybedtools import BedTool
    from maize.utils.cbook import SummaryStats

    gfffile, juncsbed = args.gff, args.bed
    tagged = "{0}.{1}.gff3".format(gfffile.rsplit(".", 1)[0], "tag_introns")

    gff3, junc = BedTool(gfffile), BedTool(juncsbed)
    ab = gff3.intersect(junc, wao=True, f=1.0, s=True)
    abfh = must_open(ab.fn)

    seen = set()
    scov = AutoVivification()

    fh, tmpgff = mkstemp(suffix=".gff3")
    fw = must_open(tmpgff, "w")
    for line in abfh:
        args = line.strip().split("\t")
        g = GffLine("\t".join(str(x) for x in args[:9]))
        if g.type == "intron" and args[10] != -1:
            ispan, jspan = g.span, int(args[11]) - int(args[10])
            if ispan == jspan:
                g.set_attr("ID", args[12], update=True)
                g.score = int(args[13])

                pparts = g.get_attr("Parent").split(".")
                locus, iso = pparts[0], ".".join(pparts[1:])
                seen.add(iso)
                if not scov[locus][iso]:
                    scov[locus][iso] = []
                scov[locus][iso].append(g.score)
            else:
                continue
        print >> fw, g
    fw.close()

    format([tmpgff, "--unique", "-o", tagged])
    os.unlink(tmpgff)

    isos = sorted(list(seen))
    fw = must_open(args.outfile, "w")
    h1, h2, stats = ["#"], ["#locus"], ["N", "mean", "median", "min", "max"]
    for iso in isos:
        h1.extend([str(iso)] + [""] * (len(stats) - 1))
        h2.extend(stats)
    print >> fw, "\t".join(str(x) for x in h1)
    print >> fw, "\t".join(str(x) for x in h2)
    for locus in scov.keys():
        out = [locus]
        for iso in isos:
            if iso in scov[locus].keys():
                juncs = scov[locus][iso]
                jstats = SummaryStats(juncs, dtype="int")
                out.extend([jstats.size, jstats.mean, jstats.median, \
                        jstats.min, jstats.max])
            else:
                out.extend(["-"] * len(stats))
        print >> fw, "\t".join(str(x) for x in out)
    fw.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'gff utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("summary", help = "print summary stats for features of different types",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input gtf file')
    sp1.add_argument("--isoform", default=False, action="store_true",
                   help="Find longest isoform of each id")
    sp1.add_argument("--ids", help="Only include features from certain IDs")
    sp1.set_defaults(func = summary)

    sp1 = sp.add_parser("filter", help = "filter the gff file based on  Identity and Coverage",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input gtf file')
    sp1.add_argument("--type", default="mRNA", help="The feature to scan for the attributes")
    g1 = sp1.add_argument_group('group1', 'Filter by identity/coverage attribute values')
    g1.add_argument("--id", default=95, type=float, help="Minimum identity")
    g1.add_argument("--coverage", default=90, type=float, help="Minimum coverage")
    g1.add_argument("--nocase", default=False, action="store_true", help="Case insensitive lookup of attribute names")
    g2 = sp1.add_argument_group('group2', 'Filter by child feature bp length')
    g2.add_argument("--child_ftype", default=None, help="Child featuretype to consider")
    g2.add_argument("--child_bp", default=None, type=int, help="Filter by total bp of children of chosen ftype")
    sp1.set_defaults(func = filter)

    sp1 = sp.add_parser('fix', help = 'fix gff fields using various options',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input GFF3 file')
    opts = 'genbank tair phytozome maize ensembl mo17 w22 ph207 phb47 hzs'.split()
    sp1.add_argument('--opt', default='ensembl', help = 'fix option')
    sp1.set_defaults(func = fix)

    sp1 = sp.add_parser('fixboundaries', help = 'fix boundaries of parent features by range chaining child features',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.add_argument("--type", default="gene", help="Feature type for which to adjust boundaries")
    sp1.add_argument("--child_ftype", default="mRNA", help="Child featuretype(s) to use for identifying boundaries")
    sp1.set_defaults(func = fixboundaries)

    sp1 = sp.add_parser('index', help = 'index gff db')
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.add_argument('fo', help = 'output GFF3 db')
    sp1.set_defaults(func = index)

    sp1 = sp.add_parser('extract', help = 'extract contig or features from gff file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.add_argument("--contigs", help="Extract features from certain contigs")
    sp1.add_argument("--names", help="Extract features with certain names")
    sp1.add_argument("--types", help="Extract features of certain feature types")
    sp1.add_argument("--children", default=0, choices=["1", "2"], help="Specify number of iterations: `1` grabs children, `2` grabs grand-children")
    sp1.add_argument("--tag", default="ID", help="Scan the tags for the names")
    sp1.add_argument("--fasta", action="store_true", help="Write FASTA if available")
    sp1.set_defaults(func = extract)

    sp1 = sp.add_parser('cluster', help = 'cluster transcripts based on shared splicing structure',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('gff', help = 'input GFF3 file')
    sp1.add_argument('outfile', help = 'output file')
    sp1.add_argument("--slop", action="store_true",
            help="allow minor variation in terminal 5'/3' UTR start/stop position")
    sp1.add_argument("--inferUTR", action="store_true",
            help="infer presence of UTRs from exon coordinates")
    sp1.set_defaults(func = cluster)

    sp1 = sp.add_parser('chain', help = 'fill in parent features by chaining children',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('gff', help = 'input GFF3 file')
    sp1.add_argument('outfile', help = 'output file')
    sp1.add_argument("--key", dest="attrib_key", default=None,
                help="Attribute to use as `key` for chaining operation")
    sp1.add_argument("--chain_ftype", default="mRNA",
                help="GFF feature type to use for chaining operation")
    sp1.add_argument("--parent_ftype", default="gene",
                help="GFF feature type to use for the chained coordinates")
    sp1.add_argument("--break", dest="break_chain", action="store_true",
                help="Break long chains which are non-contiguous")
    sp1.add_argument("--transfer_attrib", dest="attrib_list",
                help="Attributes to transfer to parent feature; accepts comma" + \
                " separated list of attribute names")
    valid_merge_op = ('sum', 'min', 'max', 'mean', 'collapse')
    sp1.add_argument("--transfer_score", dest="score_merge_op", choices=valid_merge_op,
                help="Transfer value stored in score field to parent feature." + \
                " Score is reported based on chosen operation")
    sp1.set_defaults(func = chain)

    sp1 = sp.add_parser('format', help = 'format gff file, change seqid, etc.',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('gff', help = 'input GFF3 file')
    sp1.add_argument('outfile', help = 'output file')
    valid_multiparent_ops = ["split", "merge"]
    g1 = sp1.add_argument_group('group1', 'Parameter(s) used to modify GFF attributes (9th column)')
    g1.add_argument("--name", help="Add Name attribute from two-column file")
    g1.add_argument("--note", help="Add Note attribute from two-column file")
    g1.add_argument("--add_attribute", dest="attrib_files", help="Add new attribute(s) " + \
                "from two-column file(s); attribute name comes from filename; " + \
                "accepts comma-separated list of files")
    g1.add_argument("--add_dbxref", dest="dbxref_files", help="Add new Dbxref value(s) (DBTAG:ID) " + \
                "from two-column file(s). DBTAG comes from filename, ID comes from 2nd column; " + \
                "accepts comma-separated list of files")
    g1.add_argument("--nostrict", action="store_true",
                 help="Disable strict parsing of GFF file and/or mapping file")
    g1.add_argument("--remove_attr", dest="remove_attrs", help="Specify attributes to remove; " + \
                "accepts comma-separated list of attribute names")
    g1.add_argument("--copy_id_attr_to_name", action="store_true",
                 help="Copy `ID` attribute value to `Name`, when `Name` is not defined")
    g1.add_argument("--invent_name_attr", action="store_true",
                 help="Invent `Name` attribute for 2nd level child features; " + \
                "Formatted like PARENT:FEAT_TYPE:FEAT_INDEX")
    g1.add_argument("--no_keep_attr_order", action="store_true",
                 help="Do not maintain attribute order")
    g2 = sp1.add_argument_group('group2', 'Parameter(s) used to modify content within columns 1-8')
    g2.add_argument("--seqid", help="Switch seqid from two-column file. If not" + \
                " a file, value will globally replace GFF seqid")
    g2.add_argument("--source", help="Switch GFF source from two-column file. If not" + \
                " a file, value will globally replace GFF source")
    g2.add_argument("--type", help="Switch GFF feature type from two-column file. If not" + \
                " a file, value will globally replace GFF type")
    g2.add_argument("--fixphase", action="store_true", help="Change phase 1<->2, 2<->1")
    g3 = sp1.add_argument_group('group3', 'Other parameter(s) to perform manipulations to the GFF file content')
    g3.add_argument("--unique", action="store_true", help="Make IDs unique")
    g3.add_argument("--chaindup", dest="duptype",
                 help="Chain duplicate features of a particular GFF3 `type`," + \
                      " sharing the same ID attribute")
    g3.add_argument("--multiparents", default=None, choices=valid_multiparent_ops,
                 help="Split/merge identical features (same `seqid`, `source`, `type` " + \
                 "`coord-range`, `strand`, `phase`) mapping to multiple parents")
    g3.add_argument("--remove_feats", help="Comma separated list of features to remove by type")
    g3.add_argument("--remove_feats_by_ID", help="List of features to remove by ID;" + \
                " accepts comma-separated list or list file")
    g3.add_argument("--gsac", action="store_true",
                 help="Fix GSAC GFF3 file attributes")
    g3.add_argument("--invent_protein_feat", action="store_true",
                 help="Invent a protein feature span (chain CDS feats)")
    g3.add_argument("--process_ftype", default=None,
                 help="Specify feature types to process; "
                 "accepts comma-separated list of feature types")
    g3.add_argument("--gff3", action="store_true", help="Print output in GFF3 format")
    g3.add_argument("--make_gff_store", action="store_true",
                 help="Store entire GFF file in memory during first iteration")
    sp1.set_defaults(func = format)

    sp1 = sp.add_parser('note', help = 'extract certain attribute field for each feature',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.add_argument("--type", default='gene', help="Only process certain types, multiple types allowed with comma")
    sp1.add_argument("--attribute", default="note1,note2", help="Attribute field to extract, multiple fields separated by comma")
    sp1.add_argument("--AED", type=float, help="Only extract lines with AED score <=")
    sp1.add_argument("--exoncount", action="store_true", help="Get the exon count for each mRNA feat")
    sp1.set_defaults(func = note)

    sp1 = sp.add_parser('splicecov', help = 'extract certain attribute field for each feature',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('gff', help = 'input GFF3 file')
    sp1.add_argument('bed', help = 'junctions bed')
    sp1.add_argument('outfile', help = 'output file')
    sp1.set_defaults(func = splicecov)

    sp1 = sp.add_parser('picklong', help = 'pick longest transcript')
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.set_defaults(func = pick_longest)

    sp1 = sp.add_parser('2gtf', help = 'convert gff3 to gtf format')
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.set_defaults(func = gff2gtf)

    sp1 = sp.add_parser('2tsv', help = 'convert gff3 to tsv format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.set_defaults(func = gff2tsv)

    sp1 = sp.add_parser('2bed12', help = 'convert gff3 to bed12 format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.add_argument("--parent", default="mRNA", help="Top feature type")
    sp1.add_argument("--block", default="exon", help="Feature type for regular blocks")
    sp1.add_argument("--thick", default="CDS", help="Feature type for thick blocks")
    sp1.set_defaults(func = gff2bed12)

    sp1 = sp.add_parser('2fas', help = 'extract feature (e.g. CDS) seqs and concatenate',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('gff', help = 'input GFF3 file')
    sp1.add_argument('fasta', help = 'genome fasta file')
    sp1.add_argument("--parents", dest="parents", default="mRNA",
            help="list of features to extract, use comma to separate (e.g.: 'gene,mRNA')")
    sp1.add_argument("--children", dest="children", default="CDS",
            help="list of features to extract, use comma to separate (e.g.: 'five_prime_UTR,CDS,three_prime_UTR')")
    sp1.add_argument("--feature", dest="feature",
            help="feature type to extract. e.g. `--feature=CDS` or `--feature=upstream:TSS:500`")
    valid_id_attributes = ["ID", "Name", "Parent", "Alias", "Target"]
    sp1.add_argument("--id_attribute", choices=valid_id_attributes,
            help="The attribute field to extract and use as FASTA sequence ID ")
    sp1.add_argument("--desc_attribute", default="Note",
            help="The attribute field to extract and use as FASTA sequence description")
    sp1.add_argument("--full_header", default=None, choices=["default", "tair"],
            help="Specify if full FASTA header (with seqid, coordinates and datestamp) should be generated")
    sp1.add_argument("--sep", dest="sep", default=" ", \
            help="Specify separator used to delimiter header elements")
    sp1.add_argument("--datestamp", dest="datestamp", \
            help="Specify a datestamp in the format YYYYMMDD or automatically pick `today`")
    sp1.add_argument("--conf_class", dest="conf_class", default=False, action="store_true",
            help="Specify if `conf_class` attribute should be parsed and placed in the header")
    sp1.set_defaults(func = gff2fas)

    sp1 = sp.add_parser('fromgtf', help = 'convert gtf to gff3 format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input gtf file')
    sp1.add_argument("--transcript_id", default="transcript_id", help="Field name for transcript")
    sp1.add_argument("--gene_id", default="gene_id", help="Field name for gene")
    sp1.add_argument("--augustus", default=False, action="store_true", help="Input is AUGUSTUS gtf")
    sp1.set_defaults(func = fromgtf)

    sp1 = sp.add_parser("merge", help = "merge several gff files into one",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('config', help = 'config file containing a list of gff files')
    sp1.set_defaults(func = merge)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


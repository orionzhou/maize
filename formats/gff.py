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
                    g.set_attr("ID", None)
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
    elif opt == 'as60': # strip ".1" from gene ID, no exon
        for g in gff:
            if g.type == 'gene':
                gid = g.get_attr("ID").replace(".1", '')
                g.set_attr("ID", gid)
            elif g.type.endswith('RNA'):
                gid = g.get_attr("Parent").replace(".1", '')
                g.set_attr("Parent", gid)
            elif g.type.endswith('CDS'):
                g2 = GffLine(str(g))
                g2.type = 'exon'
                g2.phase = '.'
                g2.update_attributes()
                print(g2)
            g.update_attributes()
            print(g)
    elif opt == 'evm':
        for g in gff:
            if g.type == 'gene':
                gid = g.get_attr("ID").replace("_pilon_pilon",'').replace("evm.TU.ctg", 'g')
                g.set_attr("ID", gid)
                g.set_attr("Name", None)
            elif g.type.endswith('RNA'):
                gid = g.get_attr("Parent").replace("_pilon_pilon",'').replace("evm.TU.ctg", 'g')
                rid = g.get_attr("ID").replace("_pilon_pilon",'').replace("evm.model.ctg", 'r')
                g.set_attr("ID", rid)
                g.set_attr("Parent", gid)
                g.set_attr("Name", None)
            elif g.type.endswith('CDS') or g.type.endswith("exon"):
                rid = g.get_attr("Parent").replace("_pilon_pilon",'').replace("evm.model.ctg", 'r')
                g.set_attr("Parent", rid)
                g.set_attr("ID", None)
            g.update_attributes()
            print(g)
    elif opt == 'renan': # no gene feature, 'transcipt'->'mRNA'
        for g in gff:
            if g.type.endswith('transcript'):
                mid = g.get_attr("ID")
                gid = g.get_attr("geneID")
                g.type = 'mRNA'
                g.set_attr("geneID", None)
                g.set_attr("gene_name", None)
                g.update_attributes()
                if mid.endswith('.1'):
                    g2 = GffLine(str(g))
                    g2.type = 'gene'
                    g2.set_attr("ID", gid)
                    g2.update_attributes()
                    print(g2)
                g.set_attr("Parent", gid)
            g.update_attributes()
            print(g)
    elif opt == 'nogene_noexon':
        for g in gff:
            if g.type.endswith('RNA'):
                g2 = GffLine(str(g))
                g2.type = 'gene'
                gid = g.get_attr("ID")
                mid = f'rna.{gid}'
                g.set_attr("ID", mid)
                g.set_attr("Parent", gid)
                print(g2)
            elif g.type.endswith('CDS'):
                g2 = GffLine(str(g))
                g2.type = 'exon'
                mid = "rna." + g.get_attr("Parent")
                g.set_attr("Parent", mid)
                g2.set_attr("Parent", mid)
                g2.phase = '.'
                g2.update_attributes()
                print(g2)
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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'gff utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

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

    sp1 = sp.add_parser('index', help = 'index gff db')
    sp1.add_argument('fi', help = 'input GFF3 file')
    sp1.add_argument('fo', help = 'output GFF3 db')
    sp1.set_defaults(func = index)

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

    sp1 = sp.add_parser('fromgtf', help = 'convert gtf to gff3 format',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input gtf file')
    sp1.add_argument("--transcript_id", default="transcript_id", help="Field name for transcript")
    sp1.add_argument("--gene_id", default="gene_id", help="Field name for gene")
    sp1.add_argument("--augustus", default=False, action="store_true", help="Input is AUGUSTUS gtf")
    sp1.set_defaults(func = fromgtf)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


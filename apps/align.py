#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Perform DNA-DNA alignment using BLAST, NUCMER and BLAT. Keep the interface the
same and does parallelization both in core and on grid.
"""

import os.path as op
import sys
import shutil
import logging

from maize.utils.cbook import depends
from maize.apps.base import sh, get_abs_path, which

@depends
def run_formatdb(infile=None, outfile=None, dbtype="nucl"):
    cmd = "makeblastdb"
    cmd += " -dbtype {0} -in {1}".format(dbtype, infile)
    sh(cmd)

@depends
def run_blat(infile=None, outfile=None, db="UniVec_Core",
             pctid=95, hitlen=50, cpus=16, overwrite=True):

    cmd = "pblat -threads={0}".format(cpus) if which("pblat") else "blat"
    cmd += ' {0} {1} -out=blast8 {2}'.format(db, infile, outfile)
    sh(cmd)

    blatfile = outfile
    filtered_blatfile = outfile + ".P{0}L{1}".format(pctid, hitlen)
    run_blast_filter(infile=blatfile, outfile=filtered_blatfile,
            pctid=pctid, hitlen=hitlen)
    if overwrite:
        shutil.move(filtered_blatfile, blatfile)

@depends
def run_vecscreen(infile=None, outfile=None, db="UniVec_Core",
        pctid=None, hitlen=None):
    """
    BLASTN parameters reference:
    http://www.ncbi.nlm.nih.gov/VecScreen/VecScreen_docs.html
    """
    db = get_abs_path(db)
    nin = db + ".nin"
    run_formatdb(infile=db, outfile=nin)

    cmd = "blastn"
    cmd += " -task blastn"
    cmd += " -query {0} -db {1} -out {2}".format(infile, db, outfile)
    cmd += " -penalty -5 -gapopen 4 -gapextend 4 -dust yes -soft_masking true"
    cmd += " -searchsp 1750000000000 -evalue 0.01 -outfmt 6 -num_threads 8"
    sh(cmd)

@depends
def run_megablast(infile=None, outfile=None, db=None, wordsize=None, \
        pctid=98, hitlen=100, best=None, evalue=0.01, task="megablast", cpus=16):

    assert db, "Need to specify database fasta file."

    db = get_abs_path(db)
    nin = db + ".nin"
    nin00 = db + ".00.nin"
    nin = nin00 if op.exists(nin00) else (db + ".nin")
    run_formatdb(infile=db, outfile=nin)

    cmd = "blastn"
    cmd += " -query {0} -db {1} -out {2}".format(infile, db, outfile)
    cmd += " -evalue {0} -outfmt 6 -num_threads {1}".format(evalue, cpus)
    cmd += " -task {0}".format(task)
    if wordsize:
        cmd += " -word_size {0}".format(wordsize)
    if pctid:
        cmd += " -perc_identity {0}".format(pctid)
    if best:
        cmd += " -max_target_seqs {0}".format(best)
    sh(cmd)

    if pctid and hitlen:
        blastfile = outfile
        filtered_blastfile = outfile + ".P{0}L{1}".format(pctid, hitlen)
        run_blast_filter(infile=blastfile, outfile=filtered_blastfile,
                pctid=pctid, hitlen=hitlen)
        shutil.move(filtered_blastfile, blastfile)

def run_blast_filter(infile=None, outfile=None, pctid=95, hitlen=50):
    from maize.formats.blast import filter

    logging.debug("Filter BLAST result (pctid={0}, hitlen={1})".\
            format(pctid, hitlen))
    pctidopt = "--pctid={0}".format(pctid)
    hitlenopt = "--hitlen={0}".format(hitlen)
    filter([infile, pctidopt, hitlenopt])

def main():
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'alignment utils'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("blast",
            help = 'run blastn using query against reference',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('ref_fasta', help = 'reference fasta')
    sp1.add_argument('qry_fasta', help = 'query fasta')
    task_choices = ("blastn", "blastn-short", "dc-megablast", \
                    "megablast", "vecscreen")
    sp1.add_argument("--wordsize", type=int, help="Word size")
    sp1.add_argument("--best", default=1, type=int,
            help="Only look for best N hits")
    sp1.add_argument("--task", default="megablast", choices=task_choices,
            help="Task of the blastn")
    sp1.add_argument("--pctid", type=float, default=0, help='percent identity')
    sp1.add_argument("--evalue", type=float, default=.01, help='evalue')
    sp1.add_argument("--cpus", type=int, default=1, help='cpus')
    sp1.set_defaults(func = blast)

    sp1 = sp.add_parser("blat",
            help = 'run blat using query against reference',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('ref_fasta', help = 'reference fasta')
    sp1.add_argument('qry_fasta', help = 'query fasta')
    sp1.add_argument("--pctid", type = float, default = 95)
    sp1.add_argument("--hitlen", type = int, default = 30)
    sp1.add_argument("--cpus", type=int, default=1, help='cpus')
    sp1.set_defaults(func = blast)

    sp1 = sp.add_parser("blasr",
            help = 'run blasr on a set of pacbio reads',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('ref_fasta', help = 'reference fasta')
    sp1.add_argument('fofn', help = 'fofn')
    sp1.add_argument("--cpus", type=int, default=8, help='cpus')
    sp1.set_defaults(func = blasr)

    sp1 = sp.add_parser("nucmer",
            help = 'run nucmer using query against reference',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('ref_fasta', help = 'reference fasta')
    sp1.add_argument('qry_fasta', help = 'query fasta')
    sp1.add_argument("--chunks", type=int,
                 help="Split both query and subject into chunks")
    sp1.add_argument("--cpus", type=int, default=1, help='cpus')
    sp1.add_argument("--params", default='-l 100 -c 500', help='run prarameters')
    sp1.set_defaults(func = nucmer)

    sp1 = sp.add_parser("last",
            help = 'run last using query against reference',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('query', help = 'query fasta')
    sp1.add_argument('db', help = 'subject database')
    sp1.add_argument("--dbtype", default="nucl",
                 choices=("nucl", "prot"),
                 help="Molecule type of subject database")
    sp1.add_argument("--path", help="Specify LAST path")
    sp1.add_argument("--mask", action="store_true", help="Invoke -c in lastdb")
    sp1.add_argument("--format", default="BlastTab",
                 choices=("TAB", "MAF", "BlastTab", "BlastTab+"),
                 help="Output format")
    sp1.add_argument("--minlen", default=0, type=int,
                 help="Filter alignments by how many bases match")
    sp1.add_argument("--minid", default=0, type=int, help="Minimum sequence identity")
    sp1.add_argument('-p', "--thread", type=int, default=1, help='number of threads')
    sp1.add_argument("--params", default='', help='run prarameters')
    sp1.set_defaults(func = last)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

def nucmer(args):
    """
    %prog nucmer ref.fasta query.fasta

    Run NUCMER using query against reference. Parallel implementation derived
    from: <https://github.com/fritzsedlazeck/sge_mummer>
    """
    from itertools import product
    from maize.apps.grid import MakeManager
    from maize.formats.base import split

    ref, query = args.ref_fasta, args.qry_fasta
    cpus = args.cpus
    nrefs = nqueries = args.chunks or int(cpus ** .5)
    refdir = ref.split(".")[0] + "-outdir"
    querydir = query.split(".")[0] + "-outdir"
    reflist = split([ref, refdir, str(nrefs)]).names
    querylist = split([query, querydir, str(nqueries)]).names

    mm = MakeManager()
    for i, (r, q) in enumerate(product(reflist, querylist)):
        pf = "{0:04d}".format(i)
        cmd = "nucmer -maxmatch"
        cmd += " {0}".format(args.extra)
        cmd += " {0} {1} -p {2}".format(r, q, pf)
        deltafile = pf + ".delta"
        mm.add((r, q), deltafile, cmd)
        print(cmd)

    mm.write()

def blasr(args):
    """
    %prog blasr ref.fasta fofn

    Run blasr on a set of PacBio reads. This is based on a divide-and-conquer
    strategy described below.
    """
    from maize.apps.grid import MakeManager
    from maize.utils.iter import grouper

    reffasta, fofn = args.ref_fasta, args.fofn
    flist = sorted([x.strip() for x in open(fofn)])
    h5list = []
    mm = MakeManager()
    for i, fl in enumerate(grouper(flist, 3)):
        chunkname = "chunk{0:03d}".format(i)
        fn = chunkname + ".fofn"
        h5 = chunkname + ".cmp.h5"
        fw = open(fn, "w")
        print >> fw, "\n".join(fl)
        fw.close()

        cmd = "pbalign {0} {1} {2}".format(fn, reffasta, h5)
        cmd += " --nproc {0} --forQuiver --tmpDir .".format(args.cpus)
        mm.add((fn, reffasta), h5, cmd)
        h5list.append(h5)

    # Merge h5, sort and repack
    allh5 = "all.cmp.h5"
    tmph5 = "tmp.cmp.h5"
    cmd_merge = "cmph5tools.py merge --outFile {0}".format(allh5)
    cmd_merge += " " + " ".join(h5list)
    cmd_sort = "cmph5tools.py sort --deep {0} --tmpDir .".format(allh5)
    cmd_repack = "h5repack -f GZIP=1 {0} {1}".format(allh5, tmph5)
    cmd_repack += " && mv {0} {1}".format(tmph5, allh5)
    mm.add(h5list, allh5, [cmd_merge, cmd_sort, cmd_repack])

    # Quiver
    pf = reffasta.rsplit(".", 1)[0]
    variantsgff = pf + ".variants.gff"
    consensusfasta = pf + ".consensus.fasta"
    cmd_faidx = "samtools faidx {0}".format(reffasta)
    cmd = "quiver -j 32 {0}".format(allh5)
    cmd += " -r {0} -o {1} -o {2}".format(reffasta, variantsgff, consensusfasta)
    mm.add(allh5, consensusfasta, [cmd_faidx, cmd])

    mm.write()

def get_outfile(reffasta, queryfasta, suffix="blast"):
    q = op.basename(queryfasta).split(".")[0]
    r = op.basename(reffasta).split(".")[0]
    return ".".join((q, r, suffix))

def blat(args):
    """
    %prog blat ref.fasta query.fasta

    Calls blat and filters BLAST hits.
    """
    reffasta, queryfasta = args.ref_fasta, args.qry_fasta
    blastfile = get_outfile(reffasta, queryfasta, suffix="blat")

    run_blat(infile=queryfasta, outfile=blastfile, db=reffasta,
             pctid=args.pctid, hitlen=args.hitlen, cpus=args.cpus,
             overwrite=False)

    return blastfile

def blast(args):
    """
    %prog blast ref.fasta query.fasta

    Calls blast and then filter the BLAST hits. Default is megablast.
    """
    reffasta, queryfasta = args.ref_fasta, args.qry_fasta
    blastfile = get_outfile(reffasta, queryfasta)

    run_megablast(infile=queryfasta, outfile=blastfile, db=reffasta,
                  wordsize=args.wordsize, pctid=args.pctid, evalue=args.evalue,
                  hitlen=None, best=args.best, task=args.task, cpus=args.cpus)

    return blastfile

@depends
def run_lastdb(infile=None, outfile=None, mask=False, lastdb_bin="lastdb", dbtype="nucl"):
    outfilebase = outfile.rsplit(".", 1)[0]
    db = "-p " if dbtype == "prot" else ""
    mask = "-c " if mask else ""
    cmd = "{0} {1}{2}{3} {4}".format(lastdb_bin, db, mask, outfilebase, infile)
    sh(cmd)

def last(args, dbtype=None):
    """
    %prog database.fasta query.fasta

    Run LAST by calling LASTDB and LASTAL. LAST program available:
    <http://last.cbrc.jp>

    Works with LAST-719.
    """
    query, db = args.query, args.db
    path = args.path
    nthread = args.thread
    if not dbtype:
        dbtype = args.dbtype
    getpath = lambda x: op.join(path, x) if path else x
    lastdb_bin = getpath("lastdb")
    lastal_bin = getpath("lastal")

    u = 2 if args.mask else 0
    cmd = "{0} -u {1}".format(lastal_bin, u)
    cmd += " -P {0} -i3G".format(nthread)
    cmd += " -f {0}".format(args.format)
    cmd += " {0} {1}".format(db, query)

    minlen = args.minlen
    minid = args.minid
    extra = args.params
    assert minid != 100, "Perfect match not yet supported"
    mm = minid / (100 - minid)

    if minlen:
        extra += " -e{0}".format(minlen)
    if minid:
        extra += " -r1 -q{0} -a{0} -b{0}".format(mm)
    if extra:
        cmd += " " + extra.strip()

    sh(cmd)

if __name__ == '__main__':
    main()

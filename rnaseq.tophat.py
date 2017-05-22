#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import os.path as op
import sys
import numpy as np
import argparse
import configparser
from string import Template
from colorama import init, Fore, Back, Style
import pysam

def check_bam(fbam):
    exist_bam = 1
    try:
        bam = pysam.AlignmentFile(fbam, "rb")
    except:
        exist_bam = 0
    return exist_bam
def parse_samtools_stat(fi):
    assert op.isfile(fi), "%s not exist" % fi
    r = dict()
    for line in open(fi, "r"):
        line = line.strip("\n")
        row = line.split("\t")
        if not row[0] == 'SN':
            continue
        r[row[1].strip(":")] = row[2]
    return r
def parse_tophat_alignsum(fi):
    assert op.isfile(fi), "%s not exist" % fi
    r = dict()
    for line in open(fi, "r"):
        row = line.split(":")
        if len(row) < 2:
            continue
        key = row[0].strip(" \n")
        val = row[1].strip(" \n").split(" ")[0]
        r[key] = val
    return r
def parse_tophat_alignsum_pe(fi):
    assert op.isfile(fi), "%s not exist" % fi
    r = dict()
    for line in open(fi, "r"):
        row = line.split(":")
        if len(row) < 2:
            continue
        key = row[0].strip(" \n")
        val = row[1].strip(" \n").split(" ")[0]
        if key in r:
            key += "2"
        r[key] = val
    assert r["Input"] == r["Input2"], "non-equal pairs"
    pairs = int(r["Input"])
    pairs_a = int(r['Aligned pairs'])
    pairs_ua = pairs_a - int(r["of these"])
    orphan_l = int(r["Mapped"]) - pairs_a
    orphan_r = int(r["Mapped2"]) - pairs_a
    orphan = orphan_l + orphan_r
    unmap = pairs - pairs_a - orphan 
    return [pairs, pairs_a, orphan, unmap, pairs_ua] 

def tophat(dirw, ilist, olist, diro, species, paired, 
        db_bowtie2, tophat2, samtools, samstat, parallel,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2, fo3, fo4 = "21.1.tophat.sh", "21.2.index.sh", \
            "21.3.samtools.sh", "21.4.samstat.sh"
    d22 = "22.tophat"
    fho1, fho2, fho3, fho4 = open(fo1, "w"), open(fo2, "w"), \
            open(fo3, "w"), open(fo4, "w")
    for diro in [diro]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    min_intron_len = 60000
    if species == 'rice': min_intron_len = 20000
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        fbam = "%s/%s/accepted_hits.bam" % (diro, sid)
        exist_fbam = check_bam(fbam)
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            if not exist_fbam:
                fho1.write("%s --num-threads 24 --max-multihits 20 \
                        --min-intron-length 5 --max-intron-length %d \
                        -o %s/%s %s %s %s\n" % \
                        (tophat2, min_intron_len, d22, sid, db_bowtie2, \
                        f1p, f2p))
        else:
            fr, rc, ft, rrc = row[5:9]
            if not exist_fbam:
                fho1.write("%s --num-threads 24 --max-multihits 20 \
                        --min-intron-length 5 --max-intron-length %d \
                        -o %s/%s %s %s\n" % \
                        (tophat2, min_intron_len, d22, sid, db_bowtie2, ft))
        fho2.write("samtools index %s/%s/accepted_hits.bam\n" % (diro, sid))
        fho3.write("samtools stats %s/%s/accepted_hits.bam > \
                %s/%s/samtools.stat\n" % (diro, sid, diro, sid))
        fho4.write("%s %s/%s/accepted_hits.bam\n" % (samstat, diro, sid))

    cmds = []
    cmds.append("module load bowtie2/2.2.4")
    cmds.append("cd %s" % dirw)
    cmds.append("bash %s" % fo1)
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo2))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo3))
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fo4))
    cmd = "\n".join(cmds)
    
    temdict = {
            "queue": pbs_queue,
            "walltime": pbs_walltime,
            "ppn": pbs_ppn,
            "email": pbs_email,
            "cmds": cmd
    }
    fo = "21.pbs"
    fho = open(fo, "w")
    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    fho.write(src.substitute(temdict))
    
    init()
    print(Fore.GREEN)
    print("A job script has been generated: %s" % fo)
    print("Please check, make necessary changes , then type:")
    print(Fore.RED + "qsub %s" % fo)
    print(Style.RESET_ALL)
def tophat_check(dirw, ilist, olist, diro, paired):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + [ \
                "BamFile", "MappedPairs", "OrphanPairs", "UnmappedPairs", \
                "UniquelyMappedPairs", "UniquelyMappedOrphans", \
                "InsertSizeMean", "insertSizeStd"]) + "\n")
    else:
        fho.write("\t".join(cols + [ \
                "BamFile", "Mapped", "Unmapped", "UniquelyMapped"]) + "\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        fbam = "%s/%s/accepted_hits.bam" % (diro, sid)
        assert check_bam(fbam), "%s not exist" % fbam
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            fsum = "%s/%s/align_summary.txt" % (diro, sid)
            pairs, pairs_a, orphan, unmap, pairs_ua = \
                map(lambda x: str(x), parse_tophat_alignsum_pe(fsum))
            fho.write("\t".join(row + [fbam, pairs_a, orphan, unmap, 
                pairs_ua, "", "", ""])+"\n")
        else:
            fr, rc, ft, rrc = row[5:9]
            fsum = "%s/%s/align_summary.txt" % (diro, sid)
            r1 = parse_tophat_alignsum(fsum)
            fsta = "%s/%s/samtools.stat" % (diro, sid)
            #r2 = parse_samtools_stat(fsta)
            assert rrc == r1['Input'], r1
            #assert r1['Mapped'] == r2['sequences'], (r1, r2)
            mapped = r1['Mapped']
            unmapped = str(int(rrc) - int(mapped))
            uni = str(int(mapped) - int(r1['of these']))
            fho.write("\t".join(row + [fbam, mapped, unmapped, uni])+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Map fastq to genome using tophapt2'
    )
    parser.add_argument(
            'config', nargs = '?', default = "config.ini", \
                    help = 'config file (default: config.ini)'
    )
    parser.add_argument(
            '--check', action = "store_true", \
                    help = 'run the script in check mode (default: no)'
    )
    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    cfg = cfg['tophat']
    dirw, ilist, olist, diro = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['outdir']
    species = cfg['species']
    paired = cfg.getboolean('paired')
    db_bowtie2, tophat2, samtools, samstat, parallel = \
            cfg['db_bowtie2'], cfg['tophat2'], cfg['samtools'], \
            cfg['samstat'], cfg['parallel']
    pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_template'], cfg['pbs_queue'], cfg['pbs_walltime'], \
            cfg['pbs_ppn'], cfg['pbs_email']
    if args.check:
        tophat_check(dirw, ilist, olist, diro, paired)
        sys.exit(0)
    tophat(dirw, ilist, olist, diro, species, paired, 
            db_bowtie2, tophat2, samtools, samstat, parallel,
            pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email)


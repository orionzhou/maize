#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import re
import logging
from astropy.table import Table, Column

from maize.apps.base import AttrDict, str2bool, eprint, sh, mkdir, which
from maize.formats.base import must_open
from maize.formats.pbs import PbsJob, create_job_chain

def check_cfg_index(c, noutdir = 1, njob = 3):
    c.outdirs = c.outdir.split(",")
    assert len(c.outdirs) == noutdir, "not %s outdirs: %s" % (noutdir, c.outdir)
    
    for subdir in [c.dirw, c.temp_dir] + c.outdirs:
        if not op.isdir(subdir):
            mkdir(subdir)
    
    for fn in [c.ilist, c.vcf, c.genome, c.gff]:
        assert op.isfile(fn), "cannot read %s" % fn

    for key in 'bcftools'.split():
        fp = which(c[key])
        assert fp is not None, "not executable: %s" % c[key]
        c[key] = fp
    
    c.pbs_walltimes = c.pbs_walltime.split(",")
    c.pbs_ppns = c.pbs_ppn.split(",")
    c.pbs_queues = c.pbs_queue.split(",")
    c.pbs_mems = c.pbs_mem.split(",")
    assert njob == len(c.pbs_queues) == len(c.pbs_walltimes) == len(c.pbs_ppns) == len(c.pbs_mems), "not %d jobs: %s" % (njob, c.pbs_queue)
    c.njob = njob

    return c

def check_cfg_mapping(c, noutdir = 4, njob = 2):
    c.outdirs = c.outdir.split(",")
    assert len(c.outdirs) == noutdir, "not %s outdirs: %s" % (noutdir, c.outdir)
    
    for subdir in [c.dirw, c.temp_dir] + c.outdirs:
        if not op.isdir(subdir):
            mkdir(subdir)
    
    for fn in [c.ilist, c.vcf, c.gene_bed]:
        assert op.isfile(fn), "cannot read %s" % fn

    for key in 'samtools parallel sambamba htseq bcftools bedtools'.split():
        fp = which(c[key])
        assert fp is not None, "not executable: %s" % c[key]
        c[key] = fp

    c.paired = str2bool(c.paired)
    
    if c.mapper == 'tophat2':
        c.tophat2 = which(c.tophat2)
        assert c.tophat2 is not None, "not executable: %s" % c.tophat2
    elif c.mapper == 'hisat2':
        c.hisat2= which(c.hisat2)
        assert c.hisat2 is not None, "not executable: %s" % c.hisat2
    elif c.mapper == 'star':
        c.star= which(c.star)
        assert c.star is not None, "not executable: %s" % c.star
    else:
        logging.error("unsupported mapper: %s" % c.mapper)
        sys.exit(1)
    
    assert op.isdir(c.genomedir), "cannot access %s" % c.genomedir
    genomes = set()
    t = Table.read(c.ilist, format = 'ascii.tab')
    for i in range(len(t)):
        gts = t['genome'][i].split(",")
        for gt in gts:
            genomes.add(gt)
    genomes = sorted(list(genomes))
    logging.debug("checking %d genomes" % len(genomes))
   
    c.genomes = dict()
    for gt in genomes:
        c.genomes[gt] = dict()
        dirg = "%s/%s" % (c.genomedir, gt)
        dbpre = ''
        if c.mapper == 'tophat2':
            dbpre = "%s/21.bowtie2/db" % dirg
            assert op.isfile("%s.4.bt2" % dbpre), "no %s db-index: %s" % (c.mapper, dbpre)
        elif c.mapper == 'hisat2':
            dbpre = "%s/21.hisat2/db" % dirg
            assert op.isfile("%s.8.ht2" % dbpre), "no %s db-index: %s" % (c.mapper, dbpre)
        elif c.mapper == 'star':
            dbpre = "%s/21.star" % dirg
            assert op.isfile("%s/SA" % dbpre), "no %s db-index: %s" % (c.mapper, dbpre)
        c.genomes[gt]['db'] = dbpre
        gff = "%s/51.gff" % dirg
        assert op.isfile(gff), "no gff for %s: %s" % (gff, gt)
        c.genomes[gt]['gff'] = gff
    
    c.pbs_walltimes = c.pbs_walltime.split(",")
    c.pbs_ppns = c.pbs_ppn.split(",")
    c.pbs_queues = c.pbs_queue.split(",")
    assert njob == len(c.pbs_queues) == len(c.pbs_walltimes) == len(c.pbs_ppns), "not %d jobs: %s" % (njob, c.pbs_queue)
    c.njob = njob

    return c

def index(cfg, args):
    c = AttrDict(cfg['index'])
    c = check_cfg_index(c)
    if args.check:
        return 0
    os.chdir(c.dirw)

    jcmds = [[
        "cd %s" % c.dirw
        ],[
        "cd %s" % c.dirw
        ],[
        "cd %s" % c.dirw
    ]]
    bcfgs = [
        [dict(opt = 'bash')],
        [dict(opt = 'parallel', thread = c.pbs_ppns[1])],
        [dict(opt = 'bash')]
    ]
    
    assert c.njob == len(bcfgs) == len(jcmds), "not %d jobs" % c.njob
    jobs = []
    for i in range(c.njob):
        prefix = "%s.%d" % (c.job_prefix, i+1)
        jcfg = {
            'queue': c.pbs_queues[i],
            'ppn': c.pbs_ppns[i], 
            'walltime': c.pbs_walltimes[i],
            'mem': c.pbs_mems[i], 
            'email': c.pbs_email,
        }
        job = PbsJob.from_cfg(jcfg = jcfg, jcmds = jcmds[i], bcfgs = bcfgs[i],
                prefix = prefix, njob = len(bcfgs[i]), 
                bash = c.bash, parallel = c.parallel)
        jobs.append(job)
 
    t = Table.read(c.ilist, format = 'ascii.tab')
    nrow = len(t)
    gts = [t['genotype'][x] for x in range(nrow) if t['type'][x] == 'Inbred']
    gts = set(gts)
    logging.debug("creating pseudo-refs for %d genomes" % len(gts))
    print(" ".join(gts))
    for gt in gts:
        diro = "%s/%s" % (c.outdirs[0], gt)
        mkdir(diro)
        jobs[0].subjobs[0].add_cmd("%s consensus -f %s %s -s %s \
                -c %s/25.chain -o %s/11_genome.fas" % \
                (c.bcftools, c.genome, c.vcf, gt, diro, diro))
        #jobs[1].subjobs[0].add_cmd("genome fasta %s" % diro)
        #jobs[0].subjobs[0].add_cmd("genome blat %s" % diro)
        #jobs[0].subjobs[0].add_cmd("genome bwa %s" % diro)
        jobs[1].subjobs[0].add_cmd("genome bowtie %s" % diro)
        jobs[2].subjobs[0].add_cmd("genome hisat %s" % diro)
   
    for job in jobs:
        job.write()
    fj = "%s.sh" % c.job_prefix
    create_job_chain([job.fname for job in jobs], fj)
    logging.debug("job chain with %s jobs was created: %s" % (c.njob, fj))

def mapping(cfg, args):
    c = AttrDict(cfg['mapping'])
    c = check_cfg_mapping(c)
    if args.check:
        mapping_check(c)
        return 0
    os.chdir(c.dirw)

    jcmds = [[
        "cd %s" % c.dirw,
    ], [
        #"export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % c.dirw,
    ]]
    bcfgs = [
        [dict(opt = 'bash')], 
        [dict(opt = 'parallel', thread = c.pbs_ppns[1]),
        dict(opt = 'parallel', thread = c.pbs_ppns[1])
        ],
    ]
    
    assert c.njob == len(bcfgs) == len(jcmds), "not %d jobs" % c.njob
    jobs = []
    for i in range(c.njob):
        prefix = "%s.%d" % (c.job_prefix, i+1)
        jcfg = {
            'queue': c.pbs_queues[i],
            'ppn': c.pbs_ppns[i], 
            'walltime': c.pbs_walltimes[i],
            'email': c.pbs_email,
        }
        job = PbsJob.from_cfg(jcfg = jcfg, jcmds = jcmds[i], bcfgs = bcfgs[i],
                prefix = prefix, njob = len(bcfgs[i]), 
                bash = c.bash, parallel = c.parallel)
        jobs.append(job)
 
    t = Table.read(c.ilist, format = 'ascii.tab')
    nrow = len(t)
    for i in range(nrow):
        sid = str(t['sid'][i])
        genomes = t['genome'][i].split(",")
        #logging.debug("mapping %s to %s" % (sid, ", ".join(genomes)))
        for genome in genomes:
            dbpre = c.genomes[genome]['db']
            gff = c.genomes[genome]['gff']
            pre1= "%s/%s_%s" % (c.outdirs[0], sid, genome)
            input_str = ''
            if c.paired:
                f1p = t["TrimmedReadFile1Paired"][i]
                f1u = t["TrimmedReadFile1Unpaired"][i]
                f2p = t["TrimmedReadFile2Paired"][i]
                f2u = t["TrimmedReadFile2Unpaired"][i]
                if c.mapper == 'hisat2':
                    input_str = "-1 %s -2 %s -U %s,%s" % (f1p, f2p, f1u, f2u)
                elif c.mapper == 'tophat2':
                    input_str = "%s %s,%s,%s" % (f1p, f2p, f1u, f2u)
                elif c.mapper == 'star':
                    input_str = "%s %s" % (f1p, f2p)
            else:
                ft = t["TrimmedReadFile"][i]
                if c.mapper == 'hisat2':
                    input_str = "-U %s" % ft
                elif c.mapper == 'bowtie2' or c.mapper == 'star':
                    input_str = "%s" % ft
            fbam = "%s.bam" % pre1
            if c.mapper == 'hisat2':
                jobs[0].subjobs[0].add_cmd("%s -p %s -x %s -q %s \
                        --rg-id %s --rg SM:%s | samtools view -Sb - \
                        -o %s.raw.bam" % \
                        (c.hisat2, c.pbs_ppns[0], dbpre, input_str, \
                        sid, sid, pre1))
                jobs[0].subjobs[0].add_cmd("%s sort -t %s -m 60GB %s.raw.bam -o %s.bam" % \
                        (c.sambamba, c.pbs_ppns[0], pre1, pre1))
            elif c.mapper == 'tophat2':
                jobs[0].subjobs[0].add_cmd("mkdir -p %s" % pre1)
                jobs[0].subjobs[0].add_cmd("%s -p %s -G %s \
                        --rg-id %s --rg-sample %s %s %s -o %s" % \
                        (c.tophat2, c.pbs_ppns[0], gff, sid, sid, dbpre, input_str, pre1))
                fbam = "%s/accepted.bam" % pre1
            elif c.mapper == 'star':
                jobs[0].subjobs[0].add_cmd("%s --runThreadN %s --genomeDir %s \
                        --readFilesIn %s --readFilesCommand zcat \
                        --outFileNamePrefix %s. --outSAMtype BAM SortedByCoordinate" %\
                        (c.star, c.pbs_ppns[0], dbpre, input_str, pre1))
                fbam = "%s.Aligned.sortedByCoord.out.bam" % pre1
                jobs[0].subjobs[0].add_cmd("%s index -t %s %s" % (c.sambamba, c.pbs_ppns[0], fbam))
        
            pre2 = "%s/%s_%s" % (c.outdirs[1], sid, genome)
            jobs[1].subjobs[0].add_cmd("bam stat %s --isize %s.ins.tsv > %s.tsv" % \
                    (fbam, pre2, pre2))
        
            pre3 = "%s/%s_%s" % (c.outdirs[2], sid, genome)
            fsen = "%s.txt" % pre3
            fant = "%s.as.txt" % pre3
            #if not op.isfile(fsen) or os.stat(fsen).st_size == 0:
            sam_filter_tag = "-f 1" if c.paired else ""
            #jobs[1].subjobs[1].add_cmd("%s view %s -F 256 %s | %s -r pos -s %s \
            #            -t exon -i gene_id -m union -a 20 - %s > %s" % \
            #            (c.samtools, sam_filter_tag, fbam, \
            #            c.htseq, 'reverse', gff, fsen))
            #if not op.isfile(fant) or os.stat(fant).st_size == 0:
            #jobs[2].subjobs[2].add_cmd("%s view %s -F 256 %s | %s -r pos -s %s \
            #            -t exon -i gene_id -m union -a 20 - %s > %s" % \
            #            (c.samtools, sam_filter_tag, fbam, \
            #            c.htseq, 'yes', c.gffs[j], fant))
   
    for job in jobs:
        job.write()
    fj = "%s.sh" % c.job_prefix
    create_job_chain([job.fname for job in jobs], fj)
    logging.debug("job chain with %s jobs was created: %s" % (c.njob, fj))

def mapping_check(cfg):
    t = Table.read(ilist, format = 'ascii.tab')
    nrow = len(t)
    newcols = ''
    if c.paired:
        newcols = '''BAM Pair Pair_Map Pair_Orphan Pair_Unmap
            Pair_Map_Hq Unpair Unpair_Map Unpair_Map_Hq'''.split()
    else: 
        newcols = '''BAM Total Mapped Mapped_Hq'''.split()
    for newcol in newcols:
        t.add_column(Column(name = newcol, length = nrow, dtype = object))
    
    for i in range(nrow):
        sid = t['sid'][i]
        bam = "%s/%s.bam" % (c.outdirs[0], sid)
        assert op.isfile(bam), "%s not exist" % bam
        fs = "%s/%s.tsv" % (c.outdirs[1], sid)
        assert op.isfile(fs), "%s not exist" % fs
        if c.paired:
            t['BAM'][i] = 0#bam
            t['Pair'][i] = 0#pair
            t['Pair_Map'][i] = 0#pair_map
            t['Pair_Orphan'][i] = 0#pair_orphan
            t['Pair_Unmap'][i] = 0#pair_unmap
            t['Pair_Map_Hq'][i] = 0#pair_map_hq
            t['Unpair'][i] = 0#unpair
            t['Unpair_Map'][i] = 0#unpair_map
            t['Unpair_Map_Hq'][i] = 0#unpair_map_hq
        else:
            t['BAM'] = 0#bam
            t['Total'] = 0#unpair
            t['Mapped'] = 0#unpair_map
            t['Mapped_Hq'] = 0#unpair_map_hq
    t.write(t.olist, format='ascii.tab', overwrite=True)

def run_ase(cfg):
    import pysam
    cfg = cfg['ase']
    dirw, ilist, olist, jobpre, diro = \
            cfg['dirw'], cfg['ilist'], cfg['olist'], cfg['job_prefix'], \
            cfg['outdir']
    f_fas = cfg['genome']
    paired = cfg.getboolean('paired')
    samtools, bcftools, parallel = \
            cfg['samtools'], cfg['bcftools'], cfg['parallel']
    target_vcf, gene_bed = cfg['targetvcf'], cfg['gene_bed']
    pbs_queue, pbs_walltime, pbs_ppn, pbs_email = \
            cfg['pbs_queue'], cfg['pbs_walltime'], cfg['pbs_ppn'], cfg['pbs_email']

    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    dirj = "%s.jobs" % jobpre
    if op.isdir(dirj):
        os.system("rm -rf %s" % dirj)
    for do in [diro, dirj]:
        if not op.isdir(do): 
            os.makedirs(do)
    fj = "%s.sh" % jobpre
    fhj = open(fj, "w")
    i = 1
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        gt = row[3]
        if paired:
            fbam = row[11]
        else:
            fbam = row[5]
        pre = "%s/%s" % (diro, sid)
        cmds = [
            "mkdir %s" % pre,
            "bam2bed.py %s %s.1.bed" % (fbam, pre),
            "sort -T %s -k1,1 -k2,2n %s.1.bed > %s.2.sorted.bed" % (pre, pre, pre),
            "intersectBed -wa -wb -a %s.2.sorted.bed -b %s > %s.3.bed" % (pre, target_vcf, pre),
            "sort -T %s -k4,4 -k1,1 -k2,2n %s.3.bed > %s.4.sorted.bed" % (pre, pre, pre),
            "bed.ase.py %s.4.sorted.bed %s.5.tsv %s.6.bed" % (pre, pre, pre),
            "sort -T %s -k1,1 -k2,2n %s.6.bed > %s.7.sorted.bed" % (pre, pre, pre),
            "intersectBed -wa -wb -a %s -b %s.7.sorted.bed > %s.8.bed" % (gene_bed, pre, pre),
            "bed.ase.sum.py %s.5.tsv %s.8.bed %s.tsv" % (pre, pre, pre),
            "rm %s.[1-8].*" % pre,
            "rm -rf %s" % pre,
        ]
        fo = "%s/%03d.sh" % (dirj, i)
        fho = open(fo, "w")
        fho.write("\n".join(cmds) + "\n")
        fho.close()
        i += 1
        if not op.isfile("%s.bed" % pre) or not op.isfile("%s.tsv" % pre):
            fhj.write("bash %s\n" % fo)
    fhj.close()

    cmds = []
    cmds.append("cd %s" % dirw),
    cmds.append("%s -j %s < %s" % (parallel, pbs_ppn, fj))
    
    pbsjob = PbsJob(queue = pbs_queue, 
            ppn = pbs_ppn, 
            walltime = pbs_walltime,
            email = pbs_email,
            cmds = "\n".join(cmds)
    )
    fjob = "%s.pbs" % jobpre
    pbsjob.write(fjob)
    logging.debug("Job script '%s' has been created" % fjob)

def run_ase1(dirw, ilist, olist, diro, paired, f_fas, ref_gatk,
        gatk, samtools, parallel, temp_dir,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1, fo2 = "41.1.split.sh", "41.2.hc.sh"
    fjob_pre = "41"
    fho1, fho2 = [open(x, "w") for x in [fo1, fo2]]
    for diro in [diro]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    bams = []
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        if paired:
            fbam = row[15]
        else:
            fbam = row[9]
    #    bams.append("-I %s" % fbam)
        bams.append("I=%s" % fbam)
    jgatk = "java -jar %s" % gatk #-Xmx62g
    bamstr = " ".join(bams)
    pre1 = "%s/11" % diro
    fho1.write("$PTOOL/picard.jar MergeSamFiles %s \
            O=%s.bam CREATE_INDEX=true\n" % \
            (bamstr, pre1))
    ##fho1.write("$PTOOL/picard.jar MarkDuplicates %s \
    ##        O=%s.bam M=%s.txt CREATE_INDEX=true\n" % \
    ##        (bamstr, pre1, pre1))
    pre4 = "%s/14.splitn" % diro
    fho1.write("%s -T SplitNCigarReads \
            -R %s -I %s.bam -o %s.bam -rf ReassignOneMappingQuality \
            -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS\n" % \
            (jgatk, ref_gatk, pre1, pre4))
    fho2.write("%s -T HaplotypeCaller -nct 24 \
            -R %s -dontUseSoftClippedBases \
            --maxReadsInRegionPerSample 10000 \
            -stand_call_conf 20 -I %s.bam -o %s/21.raw.vcf\n" % \
            (jgatk, ref_gatk, pre4, diro))
    #fho1.write("%s mpileup -ugf %s %s | bcftools call -vmO v -o %s" \
    #        % (samtools, f_fas, bamstr, fv))
    #fho1.write("%s -T RealignerTargetCreator -nt %s \
    #        --filter_reads_with_N_cigar \
    #        -R %s %s -o %s/01.rt.list\n" % \
    #        (jgatk, pbs_ppn, ref_gatk, bamstr, diro))
    #fho1.write("%s -T IndelRealigner -R %s \
    #        --filter_reads_with_N_cigar \
    #        %s -targetIntervals %s/01.rt.list -o %s/02.realigned.bam\n" % \
    #        (jgatk, ref_gatk, bamstr, diro, diro))
    #fho1.write("%s -T UnifiedGenotyper -nt 8 -nct 3 -R %s \
    #        --filter_reads_with_N_cigar \
    #        %s/02.realigned.bam -o %s/11.raw.vcf\n" % \
    #        (jgatk, ref_gatk, diro, diro))
    #fho1.write("vcf2tsv.py 11.raw.vcf 11.raw.tsv\n")
    
    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    
    pbs_walltimes = pbs_walltime.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_queues = pbs_queue.split(",")
    cmds = [[
        "module load picard/2.3.0",
        "module load java/jdk1.8.0_45",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "bash %s" % fo1
    ], [
        "module load picard/2.3.0",
        "module load java/jdk1.8.0_45",
        "export _JAVA_OPTIONS='-Djava.io.tmpdir=%s'" % temp_dir,
        "cd %s" % dirw,
        "bash %s" % fo2
    ]]
    njob = len(cmds)
    assert len(pbs_walltimes) == njob, "not %d jobs" % njob
    assert len(pbs_ppns) == njob, "not %d jobs" % njob

    fjobs = ["%s.%s.pbs" % (fjob_pre, chr(97+i)) for i in range(njob)]
    for i in range(njob):
        temdict = {
                "queue": pbs_queues[i],
                "walltime": pbs_walltimes[i],
                "ppn": pbs_ppns[i],
                "email": pbs_email,
                "cmds": "\n".join(cmds[i])
        }
        fho = open(fjobs[i], "w")
        fho.write(src.substitute(temdict))

    init()
    print(Fore.GREEN)
    print("%s job scripts has been generated: %s" % (njob, ", ".join(fjobs)))
    print("Please check, make necessary changes, then type:")
    print(Fore.RED + "qsub -W depend=afterany:??? %s" % fjobs[1])
    print(Style.RESET_ALL)
def run_ase2(dirw, ilist, olist, diro, paired, f_fas, target_vcf,
        samtools, bcftools, parallel,
        pbs_template, pbs_queue, pbs_walltime, pbs_ppn, pbs_email):
    if not op.isdir(dirw): os.makedirs(dirw)
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    fo1 = "24.1.bcftools.sh"
    fjob_pre = "24"
    fho1 = open(fo1, "w")
    #fho1 = [open(x, "w") for x in [fo1]]
    for diro in [diro]:
        if not op.isdir(diro): 
            os.makedirs(diro)
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        #if sid > 'BR041':
        #    continue
        if paired:
            fbam = row[15]
        else:
            fbam = row[9]
        fho1.write("%s mpileup -Ou -f %s %s | \
                %s call -C alleles -m -T %s -O v | vcf2ase.py - %s/%s.tsv\n" % 
                (bcftools, f_fas, fbam, bcftools, target_vcf, diro, sid))
    
    assert op.isfile(pbs_template), "cannot read template: %s" % pbs_template
    fht = open(pbs_template, "r")
    src = Template(fht.read())
    
    pbs_walltimes = pbs_walltime.split(",")
    pbs_ppns = pbs_ppn.split(",")
    pbs_queues = pbs_queue.split(",")
    cmds = [[
        "cd %s" % dirw,
        "%s -j %s < %s" % (parallel, pbs_ppns[0], fo1)
    ]
    ]
    njob = len(cmds)
    assert len(pbs_walltimes) == njob, "not %d jobs" % njob
    assert len(pbs_ppns) == njob, "not %d jobs" % njob

    fjobs = ["%s.%s.pbs" % (fjob_pre, chr(97+i)) for i in range(njob)]
    for i in range(njob):
        temdict = {
                "queue": pbs_queues[i],
                "walltime": pbs_walltimes[i],
                "ppn": pbs_ppns[i],
                "email": pbs_email,
                "cmds": "\n".join(cmds[i])
        }
        fho = open(fjobs[i], "w")
        fho.write(src.substitute(temdict))

    init()
    print(Fore.GREEN)
    print("%s job scripts has been generated: %s" % (njob, ", ".join(fjobs)))
    print("Please check, make necessary changes, then type:")
    #print(Fore.RED + "qsub -W depend=afterany:??? %s" % fjobs[1])
    print(Style.RESET_ALL)
def ase_check(dirw, ilist, olist, diro, paired):
    os.chdir(dirw)
    assert op.isfile(ilist), "%s not exist" % ilist
    ary = np.genfromtxt(ilist, names = True, dtype = object, delimiter = "\t")
    cols = list(ary.dtype.names)
    fho = open(olist, "w")
    if paired:
        fho.write("\t".join(cols + ["BAM",
            "Pair", "Pair_Map", "Pair_Orphan", "Pair_Unmap", \
            "Pair_Map_Hq", "Unpair", "Unpair_Map", "Unpair_Map_Hq"])+"\n")
    else:
        fho.write("\t".join(cols + ["BAM"]) + "\n")
    for row in ary:
        row = [str(x, 'utf-8') for x in list(row)]
        sid = row[0]
        bam = "%s/%s.bam" % (diro, sid)
        assert check_bam(bam), "%s not exist" % bam
        fs = "%s/%s.sum.txt" % (diro, sid)
        rs1 = picard.parse(fs)['metrics']['contents']
        rs = { rs1[i]['CATEGORY']: rs1[i] for i in list(range(len(rs1))) }
        if paired:
            f1r, f2r, rc, f1p, f1u, f2p, f2u, rrc, rc1, rc2 = row[5:15]
            pair = rs['FIRST_OF_PAIR']['TOTAL_READS']
            pair_map = rs['FIRST_OF_PAIR']['READS_ALIGNED_IN_PAIRS']
            pair_map1 = rs['FIRST_OF_PAIR']['PF_READS_ALIGNED']
            pair_map_hq1 = rs['FIRST_OF_PAIR']['PF_HQ_ALIGNED_READS']
            pair_map2 = rs['SECOND_OF_PAIR']['PF_READS_ALIGNED']
            pair_map_hq2 = rs['SECOND_OF_PAIR']['PF_HQ_ALIGNED_READS']
            unpair = rs['UNPAIRED']['TOTAL_READS']
            unpair_map = rs['UNPAIRED']['PF_READS_ALIGNED']
            unpair_map_hq = rs['UNPAIRED']['PF_HQ_ALIGNED_READS']
            pair_orphan = pair_map1 + pair_map2 - pair_map * 2
            pair_unmap = pair - pair_map - pair_orphan
            pair_map_hq = int((pair_map_hq1+pair_map_hq2)/2)
            assert pair == int(rrc), "error 1"
            assert int(rc1)+int(rc2) == unpair, "error 2"
            stats = map(str, [pair, pair_map, pair_orphan, pair_unmap, \
                    pair_map_hq, unpair, unpair_map, unpair_map_hq])
            fho.write("\t".join(row + [bam] + list(stats)) + "\n")
        else:
            fr, rc, ft, rrc = row[5:9]
            fho.write("\t".join(row + [bam]) + "\n")

    if prog == "hisat2":
        check_hisat2(cfg)
    elif prog == "htseq":
        check_htseq(cfg)

if __name__ == "__main__":
    import argparse
    import configparser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'Illumina RNA-Seq pipeline'
    )
    parser.add_argument('--config', "--cfg", default = "config.ini", help = 'config file')
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("index",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'create genome index'
    )
    sp1.add_argument("--check", action = 'store_true', help = "run script in check mode")
    sp1.set_defaults(func = index)
    
    sp1 = sp.add_parser("mapping",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = 'mapping, counting and ASE'
    )
    sp1.add_argument("--check", action = 'store_true', help = "run script in check mode")
    sp1.set_defaults(func = mapping)

    args = parser.parse_args()
    assert op.isfile(args.config), "cannot read %s" % args.config
    cfg = configparser.ConfigParser()
    cfg._interpolation = configparser.ExtendedInterpolation()
    cfg.read(args.config)
    if args.command:
        args.func(cfg, args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()


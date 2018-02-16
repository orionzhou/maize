# maize utility libraries

Collection of Python libraries to parse bioinformatics files, or perform 
common tasks related to maize assembly, annotation, and comparative genomics.

| | |
| --- | --- |
| Authors | Peng Zhou ([orionzhou](http://github.com/orionzhou)) |
| Email   | <zhoupenggeni@gmail.com> |
| License | [BSD](http://creativecommons.org/licenses/BSD/) |

## Contents

Following modules are available as generic Bioinformatics handling methods.

- `apps`
  - GenBank entrez accession, phytozome, ensembl and SRA downloader.
  - Basic phylogenetic tree construction using PHYLIP, PhyML, or RAxML, 
    and viualization.
  - Wrapper for BLAST+, LASTZ, LAST, BWA, BOWTIE2, CLC, CDHIT, CAP3, etc.

- `formats`
    Currently supports `.agp` (goldenpath), `.bed` format, `.blast` output,
    `.coords` format (`nucmer` output), `.fasta` format, `.fastq` format,
    `.gff` format, `obo` format (ontology),
    `.psl` format (UCSC blat, GMAP, etc.), `.sam` format (read mapping), etc.

- `graphics`
  - BLAST or synteny dot plot.
  - Paint regions on set of chromosomes.
  - Macro-synteny and micro-synteny plots.

- `utils`
  - Grouper can be used as disjoint set data structure.
  - range contains common range operations, like overlap
    and chaining.
  - Miscellaneous cookbook recipes, iterators decorators,
    table utilities.

Then there are modules that contain domain-specific methods.

- `annotation`
  - Calculate gene, exon and intron statistics.

- `compara`
  - Synteny scan (de-novo) and lift over (find nearby anchors).
  - Ortholog and tandem gene duplicates finder.

## Applications

Please visit [wiki](https://github.com/orionzhou/jcvi/wiki) for
full-fledged applications. 

## Dependencies

Following are a list of third-party python packages that are used by
some routines in the library. These dependencies are *not* mandatory
since they are only used by a few modules.

- [Biopython](http://www.biopython.org)
- [numpy](http://numpy.scipy.org)
- [pysam](http://pysam.readthedocs.io/en/latest)

There are other Python modules here and there in various scripts. The
best way is to install them via `pip install` when you see `ImportError`.

## Installation

The easiest way is to install it via PyPI:

```bash
pip install maize
```

To install the development version:

```bash
pip install git+git://github.com/orionzhou/maize.git
```

Alternatively, if you want to install manually:

```bash
cd ~/code  # or any directory of your choice
git clone git://github.com/orionzhou/maize.git
export PYTHONPATH=~/code:$PYTHONPATH
```

Please replace `~/code` above with whatever you like, but it must
contain `maize`. To avoid setting `PYTHONPATH` everytime, please insert
the `export` command in your `.bashrc` or `.bash_profile`.

In addition, a few module might ask for locations of external programs,
if the extended cannot be found in your `PATH`. The external programs
that are often used are:

- [bedtools](https://github.com/arq5x/bedtools2)
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [Kent tools](http://hgdownload.cse.ucsc.edu/admin/jksrc.zip)

Most of the scripts in this package contains multiple actions. To use
the `fasta` example:

```bash
usage: python -m maize.formats.fasta [-h] {summary,fas2aln,rmdot,cleanid} ...

fasta utilities

optional arguments:
  -h, --help            show this help message and exit

available commands:
  {summary,fas2aln,rmdot,cleanid}
    summary             Report length and description for each sequence
    fas2aln             convert fasta alignment file to clustal format
    rmdot               replace periods (.) in an alignment fasta by dashes
                        (-)
    cleanid             clean sequence IDs in a fasta file
```

Then you need to use one action, you can just do:

```bash
python -m maize.formats.fasta size
```

This will tell you the options and arguments it expects.

**Feel free to check out other scripts in the package, it is not just
for FASTA.**

## Reference


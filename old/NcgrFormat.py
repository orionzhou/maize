'''
Contains information about NCGR data file format and various helper methods(?).

@author: Roman Briskine
@since: 2009-12-04
'''
DELIM = '\t';
EOL = '\n';
F_REFERENCE = 0;		# Chromosome, e.g. MtChr5
F_VARIANT = 1;		# Variant ID
F_CLASS = 2;		# Polymorphism type S, D, I
F_POS = 3;			# Position on the chromosome. For indels, 263^264 or [832..842] 
F_REF_ALLELE = 4;
F_VAR_ALLELE = 5;
F_AVG_QUAL = 6;
F_MAX_QUAL = 7
F_GENE_CTXT = 8;	# Gene context: 0, I, E, U3, U5
F_REF_AA = 9;		# Reference amino acid
F_VAR_AA = 10;		# Variant amino acid
F_BLOSUM = 11;		# Blosum score
F_NEAR_GENE = 12;	# Nearest gene
F_GENE_DIST = 13;	# Gene distance
F_GENE_START = 14;
F_GENE_STOP = 15;
# Line variables
F_LN_OFFSET = 16;
F_LN_VREADS = 1;  # Number of reads calling the variant
F_LN_COV = 2;		# Coverage
F_LN_FREQ = 3;
F_LN_UNIQ = 4;		# Unique reads
F_LN_AVG_QUAL = 5;
F_LN_MAX_QUAL = 6;
F_LN_FIELDN = 7;	# Number of fields per line (accession)
# Trailing fields
F_VAR_CALLN = 1;  # Number of calls matching variant
F_REF_CALLN = 2;  # Number of calls matching reference
F_BELOWN = 3;		# Below threshold
F_NO_COV = 4;		# No coverage
F_INFO = 5;			# Matched variant info
F_TMP1 = 6;			# Redundant field names in the header
F_TMP2 = 7;
F_TMP3 = 8;
F_TMP4 = 9;
F_TRAILN = 9;		# Total number of trailing fields in the header


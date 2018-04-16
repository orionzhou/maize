'''
Parses variants file to produce a tab-delimited data file in SimpleSNP format to be used with
libsequence. The format is described in the libsequence manual. Note that tutorial gives obsolete 
information.

@author: Roman Briskine, University of Minnesota
@since: 2009-11-12
@change: 2009-12-21 Changed to work with final NCGR format
'''

import os.path;
import re;
from parsers.NcgrFormat import *;


class SimpleSnpParser():
	def __init__(self):
		self.accessionN = 0;
		self.reference = [];			# Reference alleles
		self.positions = [];			# Segregation site positions
		self.haplotypes = [];		# Variant alleles

	
	def parse(self, fPathIn, freqTheshold = 0.7, uniqThreshold = 2, fPathOut = None):
		print("Parsing...");
		if fPathOut == None:
			fPathOut = fPathIn + ".libseq";
		self.positions = [];
		self.haplotypes = [];
		
		with open(fPathIn, 'r') as fIn:
			line = fIn.readline();
			hdr = line.split(DELIM);
			# Get the variant names from the header line
			self.accessionN = (len(hdr) - F_LN_OFFSET - F_TRAILN) / F_LN_FIELDN;
			for k in range(self.accessionN):
				accNameIdx = F_LN_OFFSET + k * F_LN_FIELDN;
				self.haplotypes.append([hdr[accNameIdx]]);
			prevPos = 0;
			line = fIn.readline();
			while line != "":
				fields = line.split(DELIM);
				if fields[F_CLASS] == "S":
					if fields[F_POS] != prevPos:
						self.reference.append(fields[F_REF_ALLELE]);
						self.positions.append(fields[F_POS]);
						for k in range(self.accessionN):
							# Accession starting column
							accIdx = F_LN_OFFSET + k * F_LN_FIELDN;
							# Accession frequency
							accFreq = float(fields[accIdx + F_LN_FREQ]);
							# Accession unique reads
							accUniq = float(fields[accIdx + F_LN_UNIQ]);
							# Coverage
							accCov = int(fields[accIdx + F_LN_COV]);
							if accFreq >= freqThreshold and accUniq >= uniqThreshold:
								self.haplotypes[k].append(fields[F_VAR_ALLELE].upper());
							elif accCov >= 2:
								self.haplotypes[k].append(fields[F_REF_ALLELE]);
							else:
								self.haplotypes[k].append('N');
					else:
						for k in range(self.accessionN):
							# Accession starting column
							accIdx = F_LN_OFFSET + k * F_LN_FIELDN;
							# Accession frequency
							accFreq = float(fields[accIdx + F_LN_FREQ]);
							# Accession unique reads
							accUniq = float(fields[accIdx + F_LN_UNIQ]);
							if accFreq >= freqThreshold and accUniq >= uniqThreshold:
								self.haplotypes[k][-1] = fields[F_VAR_ALLELE].upper();
					prevPos = fields[F_POS];
				line = fIn.readline();

		with open(fPathOut, 'w') as fOut:
			#First line has the number of accessions and the number of SNP positions
			haplotypeN = len(self.haplotypes);
			snpPosN = len(self.reference);
			fOut.write(DELIM.join(map(str, [haplotypeN, snpPosN])) + EOL);
			# Segregation site positions
			fOut.write(DELIM.join(map(str, self.positions)) + EOL);
			# Reference alleles
			fOut.write(DELIM.join(self.reference) + EOL);
			# Variant alleles
			for accession in self.haplotypes:
				fOut.write(DELIM.join(accession) + EOL);


if __name__ == '__main__':
	#fPathIn = "/Users/pysar/Documents/umn/nyl/data/variants20090629/Mtchr.1.genomic_variants.20_40_2.txt";
	fPathIn = "/Users/pysar/Documents/umn/nyl/data/variants20091219/tmp.txt";
	freqThreshold = 0.7;
	uniqThreshold = 2;
	simpleSnpP = SimpleSnpParser();
	simpleSnpP.parse(fPathIn, freqThreshold, uniqThreshold);
	
	
'''
@author: Roman Briskine, University of Minnesota
'''

import os.path;
import re;


F_VARIANT = 1;
F_CLASS = 2;
F_POS = 3;
F_REF_ALLELE = 4;
F_VAR_ALLELE = 5;
F_EXON = 9;
F_ACC_OFFSET = 13;

class PhasedHaplotypeParser():
	def __init__(self, accessionN = 3, accessionColN = 7, delim = '\t'):
		self.accessionN = accessionN;
		self.accessionColN = accessionColN;
		self.delim = delim;
		self.markers = [];
		self.haplotypes = [];
		for k in range(self.accessionN + 1):
			famId = "F%03d" % k;
			self.haplotypes.append([famId]);
		self.nucleotides = { "A":1, "C":2, "G":3, "T":4 };
		

	
	def parse(self, fPathIn, freqTheshold, fPathPhased = None, fPathMarker = None):
		print("Parsing...");
		if fPathPhased == None:
			fPathPhased = fPathIn + ".haps";
		if fPathMarker == None:
			fPathMarker = fPathIn + ".info";
		
		with open(fPathIn, 'r') as fIn:
			line = fIn.readline();
			hdr = line.split(self.delim);
			self.haplotypes[0].append("REF");
			for k in range(self.accessionN):
				accNameIdx = F_ACC_OFFSET + k * self.accessionColN;
				self.haplotypes[k + 1].append(hdr[accNameIdx]);
			prevPos = 0;
			line = fIn.readline();
			while line != "":
				fields = line.split(self.delim);
				if fields[F_CLASS] == "S" and fields[F_EXON] != '' and fields[F_REF_ALLELE] in self.nucleotides:
					if fields[F_POS] != prevPos:
						self.markers.append([fields[F_VARIANT] + ":" + fields[F_EXON], fields[F_POS]]);
						nId = self.nucleotides[fields[F_REF_ALLELE]];
						self.haplotypes[0].append(nId);
						for k in range(self.accessionN):
							freqIdx = F_ACC_OFFSET + k * self.accessionColN + 3;
							if float(fields[freqIdx]) > freqThreshold:
								nId = self.nucleotides[fields[F_VAR_ALLELE].upper()];
								self.haplotypes[k + 1].append(nId);
							else:
								nId = self.nucleotides[fields[F_REF_ALLELE]];
								self.haplotypes[k + 1].append(nId);
#					else:
#						for k in range(self.accessionN):
#							freqIdx = F_ACC_OFFSET + k * self.accessionColN + 3;
#							if float(fields[freqIdx]) > freqThreshold:
#								self.haplotypes[k + 1][-1] = self.nucleotides[fields[F_VAR_ALLELE].upper()];
					prevPos = fields[F_POS];
				line = fIn.readline();

		with open(fPathMarker, 'w') as fMarker:
			for marker in self.markers:
				fMarker.write(self.delim.join(marker));
				fMarker.write('\n');
		with open(fPathPhased, 'w') as fPhased:
			for accession in self.haplotypes:
				fPhased.write(self.delim.join(map(str, accession)) + '\n');
				fPhased.write(self.delim.join(map(str, accession)) + '\n');


if __name__ == '__main__':
	fPathIn = "variant_table.10_30.txt";
	freqThreshold = 0.85;
	phParser = PhasedHaplotypeParser();
	phParser.parse(fPathIn, freqThreshold);
	
	
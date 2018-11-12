# -*- coding: utf-8 -*-
#migration
import sys;
import os;
import re;
import types;
import simuOpt;
import simuPOP;
import simuPOP.utils;
import MySQLdb;
from Bio.Seq import Seq;
from Bio.SeqRecord import SeqRecord;
from Bio import SeqIO;

dirW = os.environ["work"];
dirI = os.path.join(dirW, "Scripts/in");
dirO = os.path.join(dirW, "Scripts/out");

options = [
		{'arg': 'r:',
		 'longarg': 'n_pop=',
		 'default': 6,
		 'useDefault': True,
		 'label': 'Sub-populations',
		 'allowedTypes': [types.IntType, types.LongType],
		 'description': '',
		 'validate': simuOpt.valueBetween(0,7)
		},
		{'longarg': 'init_freq=',
		 'default': [0.2]*3 + [0.8]*3,
		 'label': 'Initial allele-freq for Sub-pop',
		 'allowedTypes': [types.ListType, types.TupleType],
		 'description': '',
		 'validate': simuOpt.valueListOf(simuOpt.valueBetween(0,1))
		},
		{'longarg': 'n_ind=',
		 'default': 100,
		 'label': 'Sub-pop Size',
		 'allowedTypes': [types.IntType, types.LongType],
		 'description': '',
		 'validate': simuOpt.valueGT(10)
		},
		{'longarg': 'n_gen=',
		 'default': 100,
		 'useDefault': True,
		 'label': 'Generations',
		 'allowedTypes': [types.IntType, types.LongType],
		 'description': '',
		 'validate': simuOpt.valueGT(0)
		},
		{'longarg': 'n_chrom=',
		 'default': 1,
		 'label': 'Chromosomes',
		 'allowedTypes': [types.IntType, types.LongType],
		 'description': '',
		 'validate': simuOpt.valueGT(0)
		},
		{'longarg': 'n_loci=',
		 'default': [1],
		 'label': 'Loci',
		 'allowedTypes': [types.ListType, types.TupleType],
		 'description': '',
		 'validate': simuOpt.valueListOf(simuOpt.valueGT(0))
		},
		{'longarg': 'n_ploidy=',
		 'default': 2,
		 'label': 'ploidy',
		 'allowedTypes': [types.IntType, types.LongType],
		 'description': '',
		 'validate': simuOpt.valueGT(0)
		},
		{'longarg': 'allele=',
		 'default': ["A", "T"],
		 'label': 'Allele',
		 'allowedTypes': [types.ListType, types.TupleType],
		 'description': '',
		 #'validate': simuOpt.valueGT(0)
		},
		#migration
		{'longarg': 'mig=',
		 'default': 0,
		 'label': 'migration model',
		 'allowedTypes': [types.IntType, types.LongType],
		 'description': '''0: No migration
		 |1: MigrIslandRates(r, n)
		 |2: MigrHierarchicalIslandRates([r11, r12], r2, [n1, n2])
		 |3: MigrSteppingStoneRates(r, n, circular=False)''',
		 'chooseOneOf': range(4),
		 'validate': simuOpt.valueOneOf(range(4))
	}
];
simuOpt.setOptions(optimized=True, alleleType='long', quiet=True)

pars = simuOpt.simuParam(options, 'A demo simulation');
cfgPath = os.path.join(dirI, 'sample.cfg');
if os.path.exists(cfgPath):
	pars.loadConfig(cfgPath);
if not pars.getParam():
	sys.exit(1);
pars.saveConfig(os.path.join(dirI, 'sample.cfg'));

n_gen,  n_pop, n_ind, n_loci, n_chrom, n_ploidy, n_allele, init_freq = \
	pars.n_gen, pars.n_pop, pars.n_ind, pars.n_loci, pars.n_chrom, pars.n_ploidy, \
	len(pars.allele), pars.init_freq;
init_freq_2 = [];
for freq in init_freq:
	init_freq_2.append([freq, 1-freq]);

pop = simuPOP.population( size=[n_ind]*n_pop, ploidy=n_ploidy, 
	loci=[n_loci], lociPos=[15, 25],
	chromNames=['chr3', 'chr5'], chromTypes=[simuPOP.Autosome], 
	alleleNames = ["A", "T"],
	infoFields='migrate_to'
	);

simu = simuPOP.simulator(pop, simuPOP.randomMating());

myInitOps, myPreOps, myPostOps = [], [], [];
#iniate sex & frequency
myInitOps += [
	simuPOP.initSex(sex=[simuPOP.Male, simuPOP.Female]),
	simuPOP.initByFreq(loci=[0], alleleFreq=init_freq_2)
];
#parameter transfer
myInitOps += [
	simuPOP.pyExec("traj_gen=[0] * n_pop"),
	simuPOP.pyExec("traj_pop=range(n_pop)"),
	simuPOP.pyExec("traj_af=init_freq"),
	simuPOP.pyExec("traj_popsize=[n_ind] * n_pop")
];
myPostOps += [
	simuPOP.stat(alleleFreq=[0], subPops=range(n_pop), vars=['alleleFreq_sp']),
	simuPOP.pyExec("traj_gen += [gen+1] * n_pop"),
	simuPOP.pyExec("traj_pop += range(n_pop)"),
	simuPOP.pyExec("traj_af += [subPop[x]['alleleFreq'][0][0] for x in range(n_pop)]"),
	simuPOP.stat(popSize=True),
	simuPOP.pyExec('traj_popsize += subPopSize')
#	simuPOP.pyEval('subPopSize'),
#	simuPOP.pyOutput("\n")
];
#iniate migration model
migMod = pars.mig;
if migMod == 1:
	myMigRate = simuPOP.utils.MigrIslandRates(0.1, n_pop);
elif migMod == 2:
	myMigRate = simuPOP.utils.MigrHierarchicalIslandRates([0.2,0.2], 0.01, [n_pop/2,n_pop/2]);
elif migMod == 3:
	myMigRate = simuPOP.utils.MigrSteppingStoneRates(0.1, n_pop, circular=True);
if migMod != 0:
	print myMigRate;
	myPreOps += [
	    #simuPOP.setAncestralDepth(2, at=[-2]),
		simuPOP.migrator(
			rate = myMigRate,
			mode = simuPOP.ByProportion,
			subPops = range(n_pop),
			toSubPops = range(n_pop)
		)
	];

simu.evolve(
	initOps = myInitOps,
	preOps = myPreOps,
	postOps = myPostOps,
	gen = n_gen
);

fPathOut = os.path.join(dirO, "out4R.dat");
fout = open(fPathOut, 'w');

gen_lst, pop_lst, af_lst, ps_lst = [], [], [], [];
for gen in simu.dvars(0).traj_gen:
	gen_lst.append(gen);
for pop in simu.dvars(0).traj_pop:
	pop_lst.append(pop);
for af in simu.dvars(0).traj_af:
	af_lst.append(af);
for ps in simu.dvars(0).traj_popsize:
	ps_lst.append(ps);

fout.write("\t".join(["gen", "pop", "af", "popsize"]) + "\n");
for i in range(len(gen_lst)):
	fout.write("%d\t%d\t%.4f\t%d\n" % (gen_lst[i], pop_lst[i], af_lst[i], ps_lst[i]));
fout.close();

#myplot = plot.myplot();
#myplot.plot(dataf);

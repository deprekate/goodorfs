#!/usr/bin/env python3
import os
import sys
import getopt

sys.path.pop(0)

from goodorfs import file_handling
from goodorfs.features import Features
from goodorfs.functions import rev_comp

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition
from sklearn.decomposition import PCA

def mad(data, axis=None):
        return np.median(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'mad', mad)
#--------------------------------------------------------------------------------------------------#
#                               ARGUMENTS                                                          #
#--------------------------------------------------------------------------------------------------#

args = file_handling.get_args()

#--------------------------------------------------------------------------------------------------#
#                               FILE INPUT                                                         #
#--------------------------------------------------------------------------------------------------#

base_trans = str.maketrans('SBVDEFHIJKLMNOPQRUWXYZ','GGGAAAAAAAAAAAAAAAAAAA')
contigs = file_handling.read_fasta(args.infile, base_trans)
if not contigs:
	sys.stdout.write("Error: no sequences found in infile\n")
	sys.exit()


#--------------------------------------------------------------------------------------------------#
#                               MAIN ROUTINE                                                       #
#--------------------------------------------------------------------------------------------------#
for id, dna in contigs.items():
	contig_features = Features(**vars(args))
	#-------------------------------Find the ORFs--------------------------------------------------#

	contig_features.parse_contig(id, dna)


	#-------------------------------Create the ORFs----------------------------------------------#
	n_unique_orfs = 0
	X = []
	Y = []
	for orfs in contig_features.iter_orfs('half'):
		for orf in orfs:
			entropy = orf.amino_acid_frequencies()
			#entropy = orf.amino_acid_entropies()
			point = []
			for aa in list('ARNDCEQGHILKMFPSTWYV#+*'):
				point.append(entropy[aa])
			point.append( orf.length() / 3)
			X.append(point)
			Y.append(orf)
		n_unique_orfs += 1

	X = StandardScaler().fit_transform(X)

	n_clust = args.kclusters #3 if n_unique_orfs<args.cutoff else 4


	#-------------------------------Cluster the ORFs----------------------------------------------#
	best = lambda: None
	best.inertia_ = float('+Inf')
	for _ in range(1000):
		model = KMeans(n_clusters=n_clust).fit(X)
		if model.inertia_ < best.inertia_:
			best = model
	labels = best.predict(X)
	
	#print("len", len(X))
	#print("mad", [np.mean(np.mad(X[labels==i,:-1], axis=0)) for i in range(n_clust)])
	#print("inertia = ", best.inertia_)

	#-------------------------------Pick the Best Cluster------------------------------------------#
	cluster = lambda : None
	cluster.mad = float('+Inf')
	cluster.minima = None
	for i in range(n_clust):
		mad = np.sum(np.mad(X[labels==i,:-1], axis=0)) 
		unique_orfs = len({orf.stop:True for orf in np.array(Y)[labels==i]})
		if mad < cluster.mad:
			if unique_orfs / n_unique_orfs > 0.05 : 
				cluster.mad = mad
				cluster.idx = i	
			else:
				cluster.minima = i
	
	#-------------------------------Output the Best Cluster----------------------------------------#
	if args.outtype != 'fna' and args.no_header:
		args.outfile.write("Sequence file = %s\n" % contig_features.infile) 
		args.outfile.write("Number of orfs = %i\n" % len(contig_features.cds)) 

	i = 1
	seen = dict()
	for label, orf in zip(labels, Y):
		if label == cluster.idx or label == cluster.minima:
			if orf.stop not in seen:
				if args.outtype == 'fna':
					args.outfile.write(">%s_orf%i [START=%s] [STOP=%s]\n" % (id, i, orf.begin(), orf.end()) )
					args.outfile.write(orf.dna)
					args.outfile.write('\n')
				elif args.outtype == 'edp':
					pass
				else:
					args.outfile.write(     str(i).rjust(5, '0'))
					args.outfile.write(orf.begin().rjust(8, ' '))
					args.outfile.write(  orf.end().rjust(8, ' '))
					args.outfile.write(  orf.frm().rjust(4, ' '))
					args.outfile.write('   0.000')
					args.outfile.write('\n')
				i += 1
				seen[orf.stop] = True


#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#


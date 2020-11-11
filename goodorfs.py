#!/usr/bin/env python
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
from goodorfs_modules import file_handling
from goodorfs_modules.features import Features
from goodorfs_modules.functions import rev_comp

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse

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

	sys.stderr.write("Sequence file = %s\n" % contig_features.infile) 
	sys.stderr.write("Number of orfs = %i\n" % len(contig_features.cds)) 
	sys.stderr.flush()
	#-------------------------------Classify the ORFs----------------------------------------------#
	X = []
	Y = []
	for orfs in contig_features.iter_orfs('half'):
		for orf in orfs:
			counts = orf.amino_acid_entropies()
			point = []
			for aa in list('ARNDCEQGHILKMFPSTWYV#+*'):
				point.append(counts[aa])
			point.append( orf.length() / 3)
			X.append(point)
			Y.append(orf)

	X = StandardScaler().fit_transform(X)

	n_clust = 2 if len(contig_features.cds)<20 else (3 if len(contig_features.cds)<args.xvar else 4)


	best = lambda: None
	best.inertia_ = float('+Inf')
	for _ in range(1000):
		model = KMeans(n_clusters=n_clust).fit(X)
		if model.inertia_ < best.inertia_:
			best = model
	labels = best.predict(X)
	#model = KMeans(n_clusters=n_clust, init='k-means++', n_init=1000, tol=1e-10, max_iter=300, algorithm='auto').fit(X)
	#labels = model.predict(X)
	#print("Inertia = ", best.inertia_)

	v = float('+Inf')
	idx = None
	loc_min = None
	for i in range(n_clust):
		stops = dict()
		variance = np.sum(np.mad(X[labels==i,], axis=0)) 
		for row in np.array(Y)[labels==i]:
			stops[row.stop] = 1
		if (variance < v):
			if (len(stops)/len(contig_features.cds) < 0.05): 
				loc_min = i
			else:
				v = variance
				idx = i	

	i = 1
	seen = dict()
	for label, orf in zip(labels, Y):
		if label == idx or label == loc_min:
			if orf.stop not in seen:
				sys.stdout.write(     str(i).rjust(5, '0'))
				sys.stdout.write(orf.begin().rjust(8, ' '))
				sys.stdout.write(  orf.end().rjust(8, ' '))
				sys.stdout.write(  orf.frm().rjust(4, ' '))
				sys.stdout.write('   0.000')
				sys.stdout.write('\n')
				#print(">%s_orf%i [START=%s] [STOP=%s]" % (id, i, orf.begin(), orf.end()) )
				#print(orf.dna)
				i += 1
			seen[orf.stop] = True

	exit()

	if contig_features.xvar % 10 != 0:
		#X = X[:,:-1]
		pca = decomposition.PCA(n_components=2).fit(X)
		x,y = pca.fit_transform(X).T
		xvector = pca.components_[0] # see 'prcomp(my_data)$rotation' in R
		yvector = pca.components_[1]
		xs = 2.2*pca.transform(X)[:,0] # see 'prcomp(my_data)$x' in R
		ys = 2.2*pca.transform(X)[:,1]
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		import matplotlib.patches as mpatches
		from matplotlib.patches import Ellipse
		fig, ax = plt.subplots(figsize=(3.93,3.49), dpi=300)
		colors = {True:'#3CC9CF', False:'#F2766E'}
		markers = {k:v for k,v in zip({0,1,2,3,4,5,6}.difference({idx}), ['o','s','H','v','^','<', '>'])}
		ax.scatter(x, y, c=[colors[x] for x in T], marker='.', linewidths=0.0, alpha=0.4, zorder=5)
		ax.scatter(x[labels==idx], y[labels==idx], facecolor='none', cmap='Spectral', s=80, linewidths=0.3, alpha=0.3, marker='d', edgecolor='black', label='Predicted', zorder=10)
		if True:
			for i in range(n_clust):
				if i != idx:
					ax.scatter(x[labels==i], y[labels==i], facecolor='none', cmap='Spectral', s=80, linewidths=0.3, alpha=0.5, marker=markers.get(i, '1'), edgecolor='green', label='Predicted', zorder=10)
		if True:
			for xx,yy,ss in zip(x,y,S):
				ax.annotate(str(ss), (xx, yy), size=3, zorder=10)
				#print(ss)
		for i in range(len(xvector)):
			#arrows project features (ie columns from csv) as vectors onto PC axes
			ax.arrow(0, 0, xvector[i]*max(xs), yvector[i]*max(ys), color='black', alpha=0.2, width=0.00001, head_width=0.0025, zorder=11)
			ax.text(xvector[i]*max(xs)*1.1, yvector[i]*max(ys)*1.1, list("ARNDCEQGHILKMFPSTWYV#+*n")[i], color='black', alpha=0.5, zorder=11)


		ax.legend(prop={'size': 6}, handles=[mpatches.Patch(color=col, label=str(lab)) for lab,col in colors.items()])
		fig.savefig('png/' + contig_features.id + '.png', bbox_inches='tight')



	exit()


	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(contig_features, shortest_path)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#


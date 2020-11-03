import sys
import itertools
from math import log2
from math import log10
from math import sqrt
from decimal import Decimal

#from .kmeans import KMeans
from .orf import CDS
from .functions import *

import numpy as np

def nrm(data, axis=None):
    x = np.array(data - np.mean(data, axis))
    d = np.linalg.norm(x, axis=1)
    return np.linalg.norm(x, axis=1)
setattr(np, 'nrm', nrm)
def amd(data, axis=None):
    return np.median(np.absolute(data - np.median(data, axis)), axis)
setattr(np, 'amd', amd)
def aad(data, axis=None):
        return np.mean(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'aad', aad)
def mad(data, axis=None):
        return np.median(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'mad', mad)
def argmid(x):
    for i, item in enumerate(x):
        if item != min(x) and item != max(x):
            return i
setattr(np, 'argmid', argmid)

class Features(list):
	"""The class holding the orfs"""
	def __init__(self, n=0, **kwargs):
		self.__dict__.update(kwargs)

		self.n = n
		self.dna = None
		self.pstop = None
		self.cds = dict()
		self.feature_at = dict()

		nucs = ['T', 'C', 'A', 'G']
		codons = [a+b+c for a in nucs for b in nucs for c in nucs]
		#amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		amino_acids = 'FFLLSSSSYY#+CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
		self.translate_codon = dict(zip(codons, amino_acids))

	def seq(self, a, b):
			return self.dna[ a-1 : b ]

	def add_orf(self, start, stop, frame, seq, rbs):
		""" Adds an orf to the factory"""
		if len(seq) < self.min_orf_len: return

		o = CDS(start, stop, frame, self)

		orfs = self.cds
		if stop not in orfs:
			orfs[stop] = dict()
			orfs[stop][start] = o
		elif start not in orfs[stop]:
			orfs[stop][start] = o
		else:
			raise ValueError("orf already defined")
		self.add_feature(o)

	def add_feature(self, feature):
		self.feature_at[ feature.as_scaled_edge()[:2] ] = feature
		self.append(feature)


	def iter_features(self, type_ = None):
		for feature in self:
			if feature.type != type_:
				continue
			yield feature

	def iter_orfs(self, kind=None):
		orfs = self.cds
		if not kind:
			for stop in orfs:
				for start in orfs[stop]:
					yield orfs[stop][start]
		elif kind == 'in':
			for stop in orfs.keys():
				keylist = list(orfs[stop].keys())
				if(orfs[stop][keylist[0]].frame > 0):
					keylist.sort()
				else:
					keylist.sort(reverse=True)
				yield (orfs[stop][start] for start in keylist)
		elif kind == 'out':
			for stop in orfs.keys():
				keylist = list(orfs[stop].keys())
				if(orfs[stop][keylist[0]].frame > 0):
					keylist.sort(reverse=True)
				else:
					keylist.sort()
				yield (orfs[stop][start] for start in keylist)
		elif kind == 'half':
			for stop in orfs.keys():
				keylist = list(orfs[stop].keys())
				if(orfs[stop][keylist[0]].frame > 0):
					keylist.sort()
				else:
					keylist.sort(reverse=True)
				keylist = keylist[ : -(-len(keylist)//2)]
				yield (orfs[stop][start] for start in keylist)

	def get_orf(self, start, stop):
		orfs = self.cds
		if stop in orfs:
			if not start:
				return orfs[stop]
			elif start in orfs[stop]:
				return orfs[stop][start]
			else:
				return lambda: None #raise ValueError("orf with start codon not found")
		else:
			return lambda: None #raise ValueError(" orf with stop codon not found")

	def get_orfs(self, stop):
		orfs = self.cds
		if stop in orfs:
			return orfs[stop]
		else:
			raise ValueError(" orf with stop codon not found")

	def get_feature(self, left, right):
		return self.feature_at.get( (left, right) , None )

	def contig_length(self):
		return len(self.dna)
	
	def end(self, frame):
		return self.contig_length() - ((self.contig_length() - (frame-1))%3)

	def num_orfs(self, kind='all'):
		count = 0
		if kind == 'all':
			for orf in self.iter_orfs():
				count += 1
		else:
			count = len(self.cds)
		return count

	def classify_orfs(self):
		from sklearn.preprocessing import StandardScaler
		from sklearn.cluster import KMeans
		from sklearn.mixture import GaussianMixture
		from sklearn.cluster import AgglomerativeClustering
		from .kmeans import KMeans as KM
		from sklearn_extra.cluster import KMedoids

		coding = dict()
		with open(self.infile.replace('fna','gb')) as fp:
			for line in fp:
				if line.startswith('     CDS '):
					if 'complement' in line:
						stop,start = line.split()[1].replace('complement(', '').replace(')', '').split('..')
					else:
						start,stop = line.split()[1].split('..')
					coding[int(stop.replace('>','').replace('<',''))] = True #start
		X = []
		Y = []
		T = []
		S = []
		uni = 0
		for orfs in self.iter_orfs('half'):
			for orf in orfs:
				if False: #not orf.has_start() or not orf.has_stop():
					continue
				counts = orf.amino_acid_entropies()
				point = []
				for aa in list('ARNDCEQGHILKMFPSTWYV#+*'):
					point.append(counts[aa])
				point.append( orf.length() / 3)
				X.append(point)
				Y.append(orf)
				T.append(coding.get(int(orf.end().replace('>','').replace('<','')), False))
				S.append(orf.end())
			uni += 1

		X = StandardScaler().fit_transform(X)

		#model = GaussianMixture(n_components=3, n_init=10, covariance_type='spherical', reg_covar=0.00001).fit(X)
		n_clust = 2 if len(self.cds)<20 else (3 if len(self.cds)<self.xvar else 4)

		'''
		best = lambda: None
		best.inertia_ = float('+Inf')
		for _ in range(1000):
			model = KMeans(n_clusters=n_clust, n_init=1, tol=1e-4, max_iter=300).fit(X)
			if model.inertia_ < best.inertia_:
				best = model
		'''

		model = KMeans(n_clusters=n_clust, n_init=1000, tol=1e-4, max_iter=300, algorithm='auto').fit(X)
		#model = KMedoids(n_clusters=n_clust, init='k-medoids++').fit(X)
		print("Inertia", model.inertia_)
		labels = model.predict(X)
		#labels = model.predict(X)
		#labels = AgglomerativeClustering(n_clusters=n_clust, linkage='ward').fit_predict(X)
		#idx = np.argmin(model.covariances_)
		#X = X[:,:-3]
		print("nrm", [np.mean(np.nrm(X[labels==i,], axis=0)) for i in range(n_clust)])
		print("amd", [np.mean(np.amd(X[labels==i,], axis=0)) for i in range(n_clust)])
		print("var", [np.mean(np.var(X[labels==i,:-3], axis=0)) for i in range(n_clust)])
		print("var", [np.mean(np.var(X[labels==i,], axis=0)) for i in range(n_clust)])
		print("mad", [np.mean(np.mad(X[labels==i,:-3], axis=0)) for i in range(n_clust)])
		print("mad", [np.mean(np.mad(X[labels==i,], axis=0)) for i in range(n_clust)])

	
		v = float('+Inf')
		idx = None
		loc_min = None
		for i in range(n_clust):
			stops = dict()
			tot = 0
			#variance = np.sum(np.mad(X[labels==i, :-3], axis=0)) 
			if self.metric == 'mad':
				variance = np.sum(np.mad(X[labels==i,], axis=0)) 
			else:
				variance = np.sum(np.var(X[labels==i,], axis=0)) 
			np.set_printoptions(linewidth=np.inf) ; print(np.round(np.mad(X[labels==i,], axis=0), 2))
			for row in np.array(Y)[labels==i]:
				stops[row.stop] = 1
				tot += 1
			print(len(stops)/len(self.cds) )
			print(len(stops), len(self.cds), tot)
			if (variance < v):
				if (len(stops)/len(self.cds) < 0.05): 
					loc_min = i
				else:
					v = variance
					idx = i	
		print('idx:' , idx)
		print('loc_min:' , loc_min)
		
		'''
		with open("/home3/katelyn/longorfs/" + self.id + ".lo") as fp:
			for line in fp:
				if line.startswith("0"):
					col = line.rstrip().split()
					if int(col[3]) > 0:
						orf = self.get_orf(int(col[1]), int(col[2])-2 )
						orf.good = True
					else:
						orf = self.get_orf(int(col[1])-2, int(col[2]) )
						orf.good = True
		
		'''
		#l = [len(X[labels==i]) for i in range(n_clust)] ; print('len', l)
		for label, orf in zip(labels, Y):
			if label == idx or label == loc_min:
				orf.good = True
				#sys.stderr.write(orf.end() + "\n")
				#print(orf.end())


		if self.xvar % 10 != 0:
			#X = X[:,:-1]
			from sklearn import decomposition
			from sklearn.decomposition import PCA
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
			fig.savefig('png/' + self.id + '.png', bbox_inches='tight')

		TP = FP = TN = FN = 0
		tp = fp = tn = fn = 0
		for orfs in self.iter_orfs('in'):
			for orf in orfs:
				if orf.good:
				#if(orf.start_codon() == 'ATG'):
					if int(orf.end().replace('<' ,'').replace('>','')) in coding:
						TP += 1
						tp += orf.length() / 3
					else:
						FP += 1
						fp += orf.length() / 3
					break
				#else:
				#	if int(orf.end().replace('<' ,'').replace('>','')) in coding:
				#		#TP += 1
				#		fn += orf.length() / 3
				#	else:
				#		#FP += 1
				#		tn += orf.length() / 3
				#	break

		print('uni', uni)
		print('T:F:X', len(coding), len(self.cds)-len(coding), len(X), sep='\t')
		FN = len(coding) - TP
		TN = len(self.cds) - TP - FP - FN
		PRECIS = TP / (TP + FP)
		RECALL = TP / (TP + FN)
		print(TP, FP, TN, FN, sep='\t')
		print(PRECIS, 'PRECIS', sep='\t')
		print(RECALL, 'RECALL', sep='\t')
		F1 = 2 * PRECIS*RECALL / (PRECIS+RECALL) if (PRECIS+RECALL) else 0
		print(F1 ,'F1', sep='\t')

		exit()

		print(tp, fp, tn, fn, sep='\t')
		precis = tp / (tp + fp)
		recall = tp / (tp + fn)
		f1 = 2 * precis*recall / (precis+recall) if (precis+recall) else 0
		print(f1 ,'f1', sep='\t')
		exit()


	def parse_contig(self, id, dna):
		self.id = id
		self.dna = dna
	
		# The dicts that will hold the start and stop codons
		stops = {1:0, 2:0, 3:0, -1:1, -2:2, -3:3}
		starts = {1:[], 2:[], 3:[], -1:[], -2:[], -3:[]}
		if dna[0:3] not in self.start_codons:
			starts[1].append(1)
		if dna[1:4] not in self.start_codons:
			starts[2].append(2)
		if dna[2:5] not in self.start_codons:
			starts[3].append(3)

		# Reset iterator and find all the open reading frames
		states = itertools.cycle([1, 2, 3])
		for i in range(1, (len(dna)-1)):
			codon = dna[i-1:i+2]
			frame = next(states)
			if codon in self.start_codons:
				starts[frame].append(i)
			elif rev_comp(codon) in self.start_codons:
				starts[-frame].append(i+2)
			elif codon in self.stop_codons:
				stop = i+2
				for start in reversed(starts[frame]):
					seq = dna[start-1:stop]
					self.add_orf(start, stop-2, frame, seq, rbs)
		
				starts[frame] = []
				stops[frame] = stop
			elif rev_comp(codon) in self.stop_codons:
				stop = stops[-frame]
				for start in starts[-frame]:
					seq = rev_comp(dna[max(0,stop-1):start])
					self.add_orf(start-2, stop, -frame, seq, rbs)
		
				starts[-frame] = []
				stops[-frame] = i
		# Add in any fragment ORFs at the end of the genome
		for frame in [1, 2, 3]:
			for start in reversed(starts[frame]):
				stop = self.end(frame)
				seq = dna[max(0,start-1):stop]
				self.add_orf(start, stop-2, frame, seq, rbs)
			start = self.end(frame)
			if rev_comp(dna[start-3:start]) not in self.start_codons:
				starts[-frame].append(self.end(frame))	
			for start in starts[-frame]:
				stop = stops[-frame]
				seq = rev_comp(dna[max(0,stop-1):start])
				self.add_orf(start-2, stop, -frame, seq, rbs)

	def first_start(self, stop):
		if stop in self:
			list = sorted(list(self[stop].keys()))
			if(list[0] < stop):
				return list[0]
			else:
				return list[-1]



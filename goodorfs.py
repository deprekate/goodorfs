#!/usr/bin/env python
import os
import sys
import getopt

#from subprocess import Popen, PIPE, STDOUT
from goodorfs_modules import file_handling
from goodorfs_modules.features import Features
from goodorfs_modules.functions import rev_comp


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
	#-------------------------------Find the ORFs----------------------------------------------#

	contig_features.parse_contig(id, dna)

	'''
	codon_counts = dict()
	for orfs in contig_features.iter_orfs("in"):
		for orf in orfs:
			for key,value in orf.codon_counts().items():
				codon_counts[key] = codon_counts.get(key, 0) + value
			break
	a = contig_features.dna.count('a')
	t = contig_features.dna.count('t')
	g = contig_features.dna.count('g')
	c = contig_features.dna.count('c')
	for key,value in sorted(codon_counts.items(), key=lambda x: x[1], reverse=True):
		print(key, value, pcodon(contig_features.dna+rev_comp(contig_features.dna), key))
	exit()
	'''

	contig_features.score_orfs()


	'''
	end_start = dict()
	with open("/home3/katelyn/projects/PHANOTATE_DATA/" + id + "/" + id + ".starts-ends") as f:
		for line in f:
			start, end = line.rstrip().split('\t')
			end_start[end] = start

	for orfs in contig_features.iter_orfs('in'):
		for i, orf in enumerate(orfs):
			if orf.end() in end_start and orf.begin() == end_start[orf.end()]:
				print(id, orf.end(), i, len(contig_features.get_orfs(orf.stop)), sep='\t')
	exit()
	'''

	# find other features

	'''
	print('#id:\t' + id)
	for orfs in contig_features.iter_orfs('in'):
		for i, orf in enumerate(orfs):
			print(i+1, orf.begin(), orf.end(), orf.start_codon(), sep='\t', end='\t')
			for aa in list('ARNDCEQGHILKMFPSTWYV'):
				print(orf.amino_acid_count(aa), end='\t')
			print()
	exit()
	'''
	'''
	scores = list()
	for orf in contig_features.iter_orfs():
			scores.append(orf.score_rbs())

	from statistics import mean
	m = mean(scores)

	counts = {'ATG':0,'GTG':0,'TTG':0}
	for orf in contig_features.iter_orfs():
		if orf.score_rbs() > m:
			counts[orf.start_codon()] += 1

	gc = (contig_features.dna.count('G') + contig_features.dna.count('C')) / contig_features.contig_length()
	print(id, gc, counts['ATG'], counts['GTG'], counts['TTG'], sep='\t')
	exit()
	'''

	#-------------------------------Write Output ----------------------------------------------#
	file_handling.write_output(contig_features, shortest_path)

#--------------------------------------------------------------------------------------------------#
#                               END                                                                #
#--------------------------------------------------------------------------------------------------#


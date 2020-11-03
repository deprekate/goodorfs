import sys
import io
import gzip
import os.path
import itertools
import textwrap
import argparse
from argparse import RawTextHelpFormatter
from math import log



def pairwise(iterable):
	a = iter(iterable)
	return zip(a, a)

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def get_args():
	usage = 'goodorfs.py [-opt1, [-opt2, ...]] infile'
	parser = argparse.ArgumentParser(description='GOODORFS: For finding and classifying open reading frames', formatter_class=RawTextHelpFormatter, usage=usage)

	parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')

	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('-s', '--start_codons', action="store", default="ATG:0.85,GTG:0.10,TTG:0.05", dest='start_codons', help='comma separated list of start codons and frequency [ATG:0.85,GTG:0.10,TTG:0.05]')
	parser.add_argument('-e', '--stop_codons', action="store", default="TAG,TGA,TAA", dest='stop_codons', help='comma separated list of stop codons [TAG,TGA,TAA]')
	parser.add_argument('-l', '--minlen', action="store", type=int, default=90, dest='min_orf_len', help='to store a variable')

	args = parser.parse_args()

	#args.start_weight = {
	#					 'ATG':1.00, 'CAT':1.00,
	#					 'GTG':0.12, 'CAC':0.12,
	#					 'TTG':0.05, 'CAA':0.05
	#					}
	start_codons = dict()
	for codon,weight in map(lambda x: tuple(x.split(':')), args.start_codons.split(',')):
		start_codons[codon.upper()] = float(weight)
	start_codons = {k: v / m for m in (max(start_codons.values()),) for k, v in start_codons.items()}
	args.start_codons = start_codons

	#args.stop_codons  = ['TAG' ,'TGA', 'TAA']
	stop_codons = []
	for codon in args.stop_codons.split(','):
		stop_codons.append(codon.upper())
	args.stop_codons = stop_codons

	return args


def read_fasta(filepath, base_trans=str.maketrans('','')):
	contigs_dict = dict()
	name = ''
	seq = ''

	lib = gzip if filepath.endswith(".gz") else io

	with lib.open(filepath, mode="rb") as f:
		for line in f:
			if line.startswith(b'>'):
				contigs_dict[name] = seq
				name = line[1:].decode("utf-8").split()[0]
				seq = ''
			else:
				#seq += line.replace("\n", "").upper()
				seq += line[:-1].decode("utf-8").upper()
		contigs_dict[name] = seq.translate(base_trans)

	if '' in contigs_dict: del contigs_dict['']
	return contigs_dict


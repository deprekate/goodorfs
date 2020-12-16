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
	usage = 'goodorfs.py [-opt1, [-opt2, ...]] infile outfile'
	parser = argparse.ArgumentParser(description='GOODORFS: For finding and classifying open reading frames', formatter_class=RawTextHelpFormatter, usage=usage)

	parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')
	#parser.add_argument('outfile', type=is_writable_file, help='output file')
	parser.add_argument('outfile', action="store", nargs='?', default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')

	parser.add_argument('-Y', '--outtype', action="store", default="tsv", dest='outtype', help='format of the output [tsv]', choices=['tsv','edp','fna'])
	parser.add_argument('-A', '--start_codons', action="store", default="ATG,GTG,TTG", dest='start_codons', help='comma separated list of start codons and frequency [ATG:0.85,GTG:0.10,TTG:0.05]')
	parser.add_argument('-Z', '--stop_codons', action="store", default="TAG,TGA,TAA", dest='stop_codons', help='comma separated list of stop codons [TAG,TGA,TAA]')
	parser.add_argument('-g', '--minlen', action="store", type=int, default=90, dest='min_orf_len', help='to store a variable')
	parser.add_argument('-n', '--no_header', action="store_false", dest='no_header', help='to store a variable')
	parser.add_argument('-k', '--kclusters', action="store", type=int, default=3, dest='kclusters', help='to store a variable')
	parser.add_argument('-c', '--cutoff', action="store", type=int, default=450, dest='cutoff', help='to store a variable')

	args = parser.parse_args()

	start_codons = []
	for codon in args.start_codons.split(','):
		start_codons.append(codon.upper())
	args.start_codons = start_codons

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


import os
import sys
from decimal import Decimal


def pcodon(seq, codon):
	length = Decimal(len(seq))
	p = dict()
	p['A'] = seq.count('A') / length
	p['T'] = seq.count('T') / length
	p['G'] = seq.count('G') / length
	p['C'] = seq.count('C') / length
	return p[codon[0]] * p[codon[1]] * p[codon[2]]

def ave(a):
	return Decimal(sum(a)/len(a))

def rev_comp(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G',
		    'N':'N',
		    'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
		    'B':'V','V':'B','D':'H','H':'D'}
	return "".join([seq_dict[base] for base in reversed(seq)])


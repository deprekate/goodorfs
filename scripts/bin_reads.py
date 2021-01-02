import sys

reads = {0:[],0.1:[],0.2:[], 0.3:[],0.4:[],0.5:[],0.6:[],0.7:[],0.8:[],0.9:[],1:[]}

def gc_content(s):
	g = s.count('G')
	c = s.count('c')
	a = s.count('A')
	t = s.count('T')
	return round( (g+c) / (g+c+a+t) , 1)

i = 0
header = ''
with open(sys.argv[1]) as fp:
	for i, line in enumerate(fp):
		if i % 2 == 0:
			header = line
		else:
			gc = gc_content(line)
			reads[gc].append(header)
			reads[gc].append(line)
			if len(reads[gc]) >= 1000:
				print("".join(reads[gc]), end='')
				print('\x00')
				reads[gc] = []

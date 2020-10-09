#!/usr/bin/env python3

import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import math

def FASTAReader(file):
	line = file.readline()

	assert line.startswith('>'), 'File format not FASTA.'

	seq_id = line[1:].strip('\n\r')

	sequence = []
	sequences = []
	line = file.readline()
	while line:
		if line.startswith('>'):
			sequences.append((seq_id, ''.join(sequence)))
			seq_id = line.rstrip('\n\r')
			sequence = []
		else:
			sequence.append(line.strip())
		line = file.readline()

	sequences.append((seq_id, ''.join(sequence)))
	return sequences

allseqs = FASTAReader(open('allseqs.mafft', 'r'))
query = allseqs[0][1]

codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

NS={}

for each in allseqs[1:]:
	refseq = each[1]
	dN, dS = 0, 0
	i = 0
	while i < len(query):
		qcodon = query[i:i+3]
		rcodon = refseq[i:i+3]
		if qcodon == rcodon or qcodon == '---' or rcodon == '---':
			i += 3
			continue
		if qcodon not in codontable.keys() or rcodon not in codontable.keys():
			i += 3
			continue
		if codontable[rcodon] == codontable[qcodon]:
			NS.setdefault(i//3, [0, 0])
			NS[i//3][1] += 1
			i += 3
			continue
		NS.setdefault(i//3, [0, 0])
		NS[i//3][0] += 1
		i += 3
		continue

df_NS = {'Pos':[], 'dN':[], 'dS':[]}
for each in NS.keys():
	df_NS['Pos'].append(each)
	df_NS['dN'].append(NS[each][0])
	df_NS['dS'].append(NS[each][1])

df_NS = pd.DataFrame(df_NS)

df_NS['dN-dS'] = df_NS['dN']-df_NS['dS']
df_NS['Z'] = df_NS['dN-dS']/df_NS['dN-dS'].std()
df_NS['P'] = st.norm.cdf(df_NS['Z'])
df_NS['p<0.05'] = (df_NS['P'] < 0.05) | (df_NS['P'] > 0.95)
df_NS['p<0.05'] = df_NS['p<0.05'].map({True: 'red', False: 'black'})

zscore = df_NS['dN-dS'].mean()/df_NS['dN-dS'].std()*math.sqrt(len(df_NS))
p = st.norm.cdf(zscore)
print('z score: ', zscore, '\n', 'p value:', p)

df_logNS = df_NS[df_NS['dN']!=0][df_NS['dS']!=0]
df_logNS['dN/dS'] = df_logNS['dN']/df_logNS['dS']
df_logNS['log2(dN/dS)'] = df_logNS['dN/dS'].apply(math.log2)

fig, ax = plt.subplots()
ax.scatter(df_logNS['Pos'], df_logNS['log2(dN/dS)'], color = df_logNS['p<0.05'], s = 5)
ax.set_xlabel('Codon position')
ax.set_ylabel('log2(dN/dS)')

fig.savefig('dNdSratio')



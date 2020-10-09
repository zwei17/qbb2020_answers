#!/usr/bin/env python3

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
			seq_id = line[1:].rstrip('\n\r')
			sequence = []
		else:
			sequence.append(line.strip())
		line = file.readline()

	sequences.append((seq_id, ''.join(sequence)))
	return sequences

dnafile = open('allseqs.fa', 'r')
alignedpepfile = open('allseqsTranslated.mafft', 'r')

dna = FASTAReader(dnafile)
alignedpep = FASTAReader(alignedpepfile)

dnafile.close()
alignedpepfile.close()

alignedDNAs = []
for dnaseq, alignpep in zip(dna, alignedpep):
	name = dnaseq[0]
	aligneddna = ''
	i=0
	for each in alignpep[1]:
		if each == '-':
			aligneddna += '---'
			continue
		aligneddna += dnaseq[1][i:i+3]
		i += 3
	alignedDNAs.append((name, aligneddna))

alignedDNAfile = open('allseqs.mafft', 'w')

for each in alignedDNAs:
	alignedDNAfile.write('>'+each[0]+'\n')
	alignedDNAfile.write(each[1]+'\n')

alignedDNAfile.close()
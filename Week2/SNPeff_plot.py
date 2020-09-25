#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

def SNPeff_plot(fname):
	vcf = open(fname, 'r')
	vcf.seek(0, 0)
	df = {}
	for line in vcf:
		if '##' in line:
			continue
		if '#' in line:
			header = line.rstrip().rsplit()
			for each in header: 
				df.setdefault(each, [])
			continue
		content = line.rstrip().rsplit(sep = '\t')
		for i in range(len(content)):
			df[header[i]].append(content[i])

	df = pd.DataFrame(df)

	DP = []
	GQ = []
	AF = []
	ANN = {}
	for each in df[fname.rsplit('.')[0] + '.fastq']:
		GQ.append(float(each.rsplit(sep = ':')[1]))
		DP.append(int(each.rsplit(sep = ':')[2]))

	for each in df['INFO']:
		info = each.rstrip().rsplit(sep = ';')
		AF.append(float(info[3].rsplit('=')[1]))
		for content in info:
			if content[0:4] == 'ANN=':
				snpeff_ann = content.rsplit(sep = ',')
				for eachann in snpeff_ann:
					eff = eachann.rsplit(sep = '|')[1]
					ANN.setdefault(eff, 0)
					ANN[eff] += 1
				break
	df['DP'] = DP
	df['GQ'] = GQ
	df['AF'] = AF
	ANN = pd.Series(ANN)

	fig, ax = plt.subplots(figsize = (10, 10), nrows = 2, ncols = 2)
	ax[0, 0].hist(df['DP'], bins = range(0, 800, 5), log = True)
	ax[0, 0].set_title('Read depth distribution')
	ax[0, 0].set_xlabel('Coverage')
	ax[0, 0].set_ylabel('Counts')
	ax[0, 1].hist(df['GQ'], bins = range(0, 200, 5))
	ax[0, 1].set_title('Read quality distribution')
	ax[0, 1].set_xlabel('Phred value')
	ax[0, 1].set_ylabel('Counts')
	ax[1, 0].hist(df['AF'], bins = 20, log = True)
	ax[1, 0].set_title('Alleilic frequency distribution')
	ax[1, 0].set_xlabel('Alleilic frequency')
	ax[1, 0].set_ylabel('Counts')
	ax[1, 1].bar(ANN.keys(), ANN, log = True)
	ax[1, 1].set_title('SNPEFF annotation distribution')
	ax[1, 1].set_ylabel('Counts')
	ax[1, 1].set_xticklabels(ANN.keys(), rotation = 90)
	plt.savefig(fname.rsplit('.')[0] + '.png', bbox_inches='tight')

	vcf.close()

	return fig

if __name__ == '__main__':
	SNPeff_plot(sys.argv[1])

#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def Manhattan_and_QQ(fname, outdir):
	assocLinear = pd.read_csv(fname, sep = '\s+')
	assocLinear['-logP'] = -1*np.log10(assocLinear['P'])
	assocLinear['idx'] = range(len(assocLinear))
	assocLinear['P<e-5?'] = assocLinear['-logP'] > 5

	fig1, ax1 = plt.subplots(figsize = (8., 6.))
	for chrs in assocLinear['CHR'].unique():
		ax1.scatter(assocLinear['idx'][assocLinear['TEST'] == 'ADD'][assocLinear['CHR'] == chrs][assocLinear['P<e-5?'] == False], assocLinear['-logP'][assocLinear['TEST'] == 'ADD'][assocLinear['CHR'] == chrs][assocLinear['P<e-5?'] == False], marker = '.', color = '#555555' ,alpha = 0.2)
		ax1.scatter(assocLinear['idx'][assocLinear['TEST'] == 'ADD'][assocLinear['CHR'] == chrs][assocLinear['P<e-5?']], assocLinear['-logP'][assocLinear['TEST'] == 'ADD'][assocLinear['CHR'] == chrs][assocLinear['P<e-5?']], marker = '.', color = '#FF0000' ,alpha = 0.2)

	plt.xlabel("SNPs")
	plt.ylabel("-log10(p)")
	plt.title('Manhattan plot for trait'+fname.rsplit(sep = '.')[1])
	
	fig1.savefig(outdir+'/'+fname.rsplit(sep = '.')[1]+'_Manhattan')

	assocSorted = assocLinear.sort_values(by = 'P')
	assocSorted['idx'] = range(1, len(assocSorted)+1)
	assocSorted['Expected p'] = assocSorted['idx']/len(assocSorted)
	assocSorted['Expected -log(p)'] = -1*np.log10(assocSorted['Expected p'])

	fig2, ax2 = plt.subplots(figsize = (8., 8.))

	ax2.scatter(assocSorted['Expected -log(p)'], assocSorted['-logP'])
	ax2.plot([int(assocSorted['Expected -log(p)'].max()), 0], [int(assocSorted['Expected -log(p)'].max()), 0], color = "black")

	plt.xlabel("Expected -log10(p)")
	plt.ylabel("Actual -log10(p)")
	plt.title('QQ plot for trait'+fname.rsplit(sep = '.')[1])
	
	fig2.savefig(outdir+'/'+fname.rsplit(sep = '.')[1]+'QQ')
	return

if __name__ == '__main__':
	Manhattan_and_QQ(sys.argv[1], sys.argv[2])

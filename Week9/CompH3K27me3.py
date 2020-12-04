#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import pyBigWig as pbw

#Load HeK27me3 data.
bw = pbw.open('data/WT_H3K27me3.bw')

#Load gene expression and compartmentalization data.
Nf = open('Negcompgenes.bed', 'r')
Pf = open('Poscompgenes.bed', 'r')
Bgenes = np.loadtxt(Nf, str, delimiter = '\t')
Agenes = np.loadtxt(Pf, str, delimiter = '\t')
Agenes = Agenes[Agenes[:, 0]=='chr13']
Nf.close()
Pf.close()

#Extract expression data and H3K27me3 data.
Acomp, Bcomp = [], []

for gene in Agenes:
	h3k27me3 = bw.stats(gene[0], int(gene[1]), int(gene[2]), type='sum')[0]
	if h3k27me3 == None:
		continue
	Acomp.append([float(gene[4]), h3k27me3])

for gene in Bgenes:
	h3k27me3 = bw.stats(gene[0], int(gene[1]), int(gene[2]), type='sum')[0]
	if h3k27me3 == None:
		continue
	Bcomp.append([float(gene[4]), h3k27me3])

Acomp = np.array(Acomp, dtype = 'float')
Bcomp = np.array(Bcomp, dtype = 'float')
Acomp = np.log10(Acomp + 1)
Bcomp = np.log10(Bcomp + 1)

#Generate scatter plots.
fig, ax = plt.subplots(figsize = (10, 4.5), nrows = 1, ncols = 2, sharey = True)
ax[0].scatter(Acomp[:, 0], Acomp[:, 1], s = 6, c = 'red', alpha = 0.5)
ax[1].scatter(Bcomp[:, 0], Bcomp[:, 1], s = 6, c = 'blue', alpha = 0.5)

ax[0].set_xlabel('log10(fpkm+1)')
ax[1].set_xlabel('log10(fpkm+1)')
ax[0].set_ylabel('log10(H3K27me3+1)')

ax[0].set_title('A compartment genes')
ax[1].set_title('B compartment genes')
plt.show()
fig.savefig('CompExpH3K27me3')

#!/usr/bin/env python2

import hifive as hf
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.collections as collections

hic = hf.HiC('mm9HiCProject', 'r')
#Load Hi-C project.

Comp = hf.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
Comp.write_eigen_scores('hic_comp.bed')
#Calculate one eigenvector of the correlation matrix.

X = Comp.positions['chr13']/100000
#Store bin numbers in X.
Y = Comp.eigenv['chr13']
#Store compartment scores in Y.

fig, ax = plt.subplots(figsize = (7, 5))
ax.plot(X, Y, color = 'black', lw = 0.5)
ax.axhline(0, color = 'black', lw = 1)

collection = collections.BrokenBarHCollection.span_where(X[:,0], ymin = 0, ymax = float(np.max(Y)), where = Y > 0, facecolor = 'red', alpha  =0.5)
ax.add_collection(collection)
collection = collections.BrokenBarHCollection.span_where(X[:,0], ymin = float(np.min(Y)), ymax = 0, where = Y < 0, facecolor = 'blue', alpha  =0.5)
ax.add_collection(collection)

ax.set_xlabel('Chr13 (x100 kb)')
ax.set_ylabel('Compartment score')
ax.set_title('Hi-C compartmentalization')

plt.show()
fig.savefig('Chr13Compartments')
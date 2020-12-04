#!/usr/bin/env python2

import hifive as hf
import numpy as np
from matplotlib import pyplot as plt

hic = hf.HiC('mm9HiCProject', 'r')
#Load Hi-C project.

data = hic.cis_heatmap('chr13', 1000000, datatype = 'fend', arraytype = 'full', diagonalincluded = True)
#Generate actual counts and expected counts.

enrichment = data[:, :, 0] / data[:, :, 1]
#Calculate enrichment scores.

pseudoenrich = np.where(enrichment >= 0, enrichment, 1)
#Replace Nan by pseudocount 1.

hicmap = np.log10(pseudoenrich + 0.5)
#Log10 transform the data.

#Plot the heatmap.
fig, ax = plt.subplots(figsize = (8, 8))
img = ax.imshow(hicmap, cmap='Reds')
fig.colorbar(img, fraction = 0.02, pad = 0.02)

ax.xaxis.set_ticks_position('top')
ax.set_xlabel('Chr13 (x100 kb)')
ax.xaxis.set_label_position('top')
ax.set_ylabel('Chr13 (x100 kb)')
ax.set_title('Chr13 Hi-C map', pad = 50)

plt.show()
fig.savefig('Chr13HiC')

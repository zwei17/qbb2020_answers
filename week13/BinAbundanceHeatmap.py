#!/usr/bin/env python3

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

binabun = pd.read_csv('abundance_table.tab', sep = '\t')
binabun = binabun.loc[:, 'SRR492186':]

binabun = np.log10(binabun + 1)
binabun = binabun.sort_index(axis = 1)

fig, ax = plt.subplots(figsize = (10, 5))
img = ax.imshow(binabun, cmap = 'hot')
fig.colorbar(img, fraction = 0.02, pad = 0.02)

ax.set_xticks(np.arange(len(binabun.columns)) - 0.5)
ax.set_xticklabels(binabun.columns)
ax.set_yticklabels(['', 'S. haemolyticus', 'E. faecalis', 'L. citreum', 'C. avidum', 'S. epidermidis', 'S. aureus subsp.', 'A. prevotii', 'S. lugdunensis'])
ax.set_title('Bin Abundance Heatmap')
plt.xticks(rotation = 45)
plt.show()
fig.savefig('BinAbundance')
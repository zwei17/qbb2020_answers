#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np

#Load gene lists into numpy arrays.
Nf = open('Negcompgenes.bed', 'r')
Pf = open('Poscompgenes.bed', 'r')
Neg = np.loadtxt(Nf, str, delimiter = '\t')
Pos = np.loadtxt(Pf, str, delimiter = '\t')
Pos = Pos[Pos[:, 0]=='chr13']
Nf.close()
Pf.close()

#Extract fpkm data and log-transform them.
N = np.array(Neg[:, 4], dtype = 'float')
P = np.array(Pos[:, 4], dtype = 'float')
N = np.log10(N + 1)
P = np.log10(P + 1)
#Mask 0 values.
#N = np.ma.masked_equal(N, 0)
#P = np.ma.masked_equal(P, 0)

#Generate violin plots.
fig, ax = plt.subplots(figsize = (8, 6))
img = ax.violinplot([P, N], positions = [1.35, 1.65], showmedians = True, widths = 0.25)

img['cbars'].set_color('black')
img['cbars'].set_linewidth(1)
img['cmedians'].set_color('black')
img['cmaxes'].set_color('black')
img['cmins'].set_color('black')

img['bodies'][1].set_facecolor('lightskyblue')
img['bodies'][1].set_edgecolor('darkblue')
img['bodies'][1].set_alpha(0.8)
img['bodies'][0].set_facecolor('salmon')
img['bodies'][0].set_edgecolor('darkred')
img['bodies'][0].set_alpha(0.8)

ax.set_ylabel('log10(fpkm+1)')
ax.get_xaxis().set_tick_params(direction='out')
ax.set_xticks([1.35, 1.65])
ax.set_xticklabels(['A compartment (positive)', 'B compartment (negative)'])
ax.set_title('Gene expression in different compartments')

plt.show()
fig.savefig('CompGeneExp')

#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

G1Efeature = pd.read_csv('G1EFeatureCount.csv', delim_whitespace = True, names = ['Count', 'Feature'])
ER4feature = pd.read_csv('ER4FeatureCount.csv', delim_whitespace = True, names = ['Count', 'Feature'])

diffpeaks = pd.read_csv('Diffpeaks.csv', delim_whitespace = True, names = ['Count', 'Feature'])

fig, ax = plt.subplots(ncols = 2, figsize = (10, 8))
ax[0].bar([-0.2, 0.8, 1.8], G1Efeature['Count'], width = 0.4, color = '#00BBAB', alpha = 0.6, label = 'G1E')
ax[0].bar([0.2, 1.2, 2.2], ER4feature['Count'], width = 0.4, color = '#FF0000', alpha = 0.6, label = 'ER4')
ax[0].set_ylabel('Count number', size = 15)
ax[0].set_xticks([0,1,2])
ax[0].set_xticklabels(G1Efeature['Feature'])
ax[0].tick_params(labelsize = 15)
ax[0].set_title('Peak feature distribution', size = 15)
ax[0].legend(loc = 2)
ax[1].bar(['Lost peaks', 'Gained peaks'], diffpeaks['Count'][:2], width = 0.6, color = '#AAAAAA')
ax[1].tick_params(labelsize = 15)
ax[1].set_title('Differential binding peaks', size = 15)

fig.savefig('PeakPlot')
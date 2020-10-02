#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np

df_GT = pd.read_csv('BYxRM_segs_saccer3.SNPcounts.raw', sep = ' ')
pca_input = df_GT.iloc[:, 6:]

pca_input_stded = StandardScaler().fit_transform(pca_input)
pca_input_stded = np.nan_to_num(pca_input_stded, nan=0)

pca = PCA(n_components = 10)
pca_output = pca.fit_transform(pca_input_stded)

pca_output_df = pd.DataFrame(data = pca_output, columns = ['PC1', 'PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'])

fig1, ax1 = plt.subplots(figsize = (10., 6.))
ax1.scatter(x = pca_output_df['PC1'], y = pca_output_df['PC2'])
ax1.set_title('PCA analysis of yeast strains')
ax1.set_xlabel('PC1')
ax1.set_ylabel('PC2')

fig1.savefig('PCA_results')

df_freq = pca_input.T
df_freq['frequency'] = df_freq.mean(axis = 1)/2

fig2, ax2 = plt.subplots(figsize = (6., 4.))
ax2.hist(df_freq['frequency'], bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], log = True)
ax2.set_title('Allele frequency spectrum')
ax2.set_xlabel('Allele frequency')
ax2.set_ylabel('Counts')

fig2.savefig('Frequency_spectrum')








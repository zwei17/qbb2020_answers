#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

GT_snp = pd.read_csv('ChrII795756.raw', sep = ' ')
GT_snp['FIID'] = GT_snp['FID'].map(str)+'_'+GT_snp['IID'].map(str)

PT = pd.read_csv('Phenotypes.txt', sep = '\t')
PT['FIID'] = PT['FID'].map(str)+'_'+PT['IID'].map(str)

df_GP = pd.merge(GT_snp, PT, on = 'FIID')
df_GP = df_GP.dropna(subset=['Cadmium_Chloride'])

fig, ax = plt.subplots(figsize = (8., 5.))
ax.boxplot([df_GP['Cadmium_Chloride'][df_GP['._C']==2], df_GP['Cadmium_Chloride'][df_GP['._C']==1], df_GP['Cadmium_Chloride'][df_GP['._C']==0]], labels = ['CC', 'CG', 'GG'])
ax.set_ylabel('Relative colocy size on cadmium chloride plate')
ax.set_xlabel('Genotypes at Chr II: 795756')
fig.savefig('Cadmium_Chloride_ChrII795756_boxplot')
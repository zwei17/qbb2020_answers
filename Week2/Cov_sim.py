#!/usr/bin/env python3

import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import sys

def Cov_sim(cov):
	sequence = []
	reads = {}

	for i in range(1, 1000001):
		reads[i] = 0
	for read in range(int(cov*10000)):
		start = random.randint(1, 999900)
		for bp in range(100):
			sequence.append(start+bp)
			reads[start+bp] += 1

	reads = pd.Series(reads)
	data = np.arange(max(reads)+8)
	sequence = pd.Series({'reads':sequence})

	fig, ax = plt.subplots(figsize = (10, 5), ncols = 2)
	ax[0].hist(sequence, bins = 10000)
	ax[0].set_xlabel('nt')
	ax[0].set_ylabel('Read counts')
	histo = ax[1].hist(reads, bins = range(0, max(reads)+8, 1), density = True)
	ax[1].plot(data, stats.poisson.pmf(data, mu = cov), color = 'black')
	ax[1].set_xlabel('Coverage')
	ax[1].set_ylabel('Frequency')
	plt.savefig('Cov'+str(cov)+ '.png', bbox_inches='tight')
	print(histo[0][0]*cov*10000)

	return fig

if __name__ == '__main__':
	Cov_sim(sys.argv[1])
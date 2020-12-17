#!/usr/bin/env python3

import sys

def KrakenToKrona(fname):
	'''This function transform .kraken files to krona chart .txt files that are recognizable by Krona tools.'''
	f = open(fname, 'r')

	data = dict()
	for line in f:
		taxo = line.rstrip().rsplit(sep = 'root;')
		if len(taxo) != 2:
			taxo.append('')
		data.setdefault(taxo[1], 0)
		data[taxo[1]] += 1 #Calculate number of reads for each specie.

	f.close()

	orgs = list(data.keys())
	out = open(fname.rsplit('.k')[0]+'.kronachart', 'w')
	for each in orgs:
		#Krona tools require that each level is separated by tabs, so re-write the dictionary into a tab-separated chart.
		content = each.rsplit(sep = ';')
		out.write(str(data[each]))
		for entry in content:
			out.write('\t')
			out.write(entry)
		out.write('\n')
    
	out.close()
	return

if __name__ == '__main__':
	for file in sys.argv[1:]:
		KrakenToKrona(file)
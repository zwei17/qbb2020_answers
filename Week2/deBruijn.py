#!/usr/bin/env python3

import sys

def deBruijn(fname):
	seqs = open(fname,'r')
	seqs.seek(0, 0)

	reads = []
	for line in seqs: reads.append(line.rstrip())
	seqs.close()

	startk = {}
	endk = {}
	for each in reads:
		startk.setdefault(each[:3], [])
		startk[each[:3]].append(each)
		endk.setdefault(each[-3:], [])
		endk[each[-3:]].append(each)

	connect = {}
	for each in reads:
		endk.setdefault(each[:3], False)
		startk.setdefault(each[-3:],False)
		connect.setdefault(each, (endk[each[:3]], startk[each[-3:]]))

	used = []
	contigs = []
	contig = ''
	start = []
	end = []
	for each in reads:
		if not connect[each][0] or len(connect[each][0]) > 1: 
			start.append(each)
		if not connect[each][1] or len(connect[each][1]) > 1:
			end.append(each)
			continue
		if len(connect[connect[each][1][0]][0]) > 1:
			end.append(each)
	for each in start:
		contig = each
		while each not in used and each not in end:
			contig += connect[each][1][0][-2:]
			used.append(each)
			each = connect[each][1][0]
		if each not in used:
			used.append(each)
		contigs.append(contig)
	for each in reads:
		if each not in used:
			contigs.append(each)
	print(contigs)
	contigfile = open('Contigs_'+fname, 'w')
	for i in range(len(contigs)):
		contigfile.write('>Contig #'+str(i)+':\n')
		contigfile.write(contigs[i]+'\n')
	contigfile.close()

	return contigs


if __name__ == '__main__':
	deBruijn(sys.argv[1])
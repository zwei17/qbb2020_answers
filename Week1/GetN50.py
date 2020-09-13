#!/usr/bin/env python3

import sys

def GetN50(Contigs, Ref_len):
	contigs = open(Contigs, 'r')
	contigs.seek(0, 0)
	covered = 0
	for line in contigs:
		line = line.rstrip().split()
		N50 = int(line[1])
		covered += N50
		if covered >= (int(Ref_len)/2.0):
			break
	contigs.close()
	return N50


if __name__ == '__main__':
	print(GetN50(sys.argv[1], sys.argv[2]))

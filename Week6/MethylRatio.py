#!/usr/bin/env python3

#Load genomic coordinates
Coord = open('mm10_refseq_genes_chr6_50M_60M.bed', 'r')
Coord.seek(0, 0)

#Load .bedGraph files
E40ICM = open('SRR3083926_1.chr6_bismark_bt2_pe.bedGraph', 'r')
E55epi = open('SRR3083929_1.chr6_bismark_bt2_pe.bedGraph', 'r')

#Create output file
output = open('MethylRatio.bed', 'w')

#Calculate methylation change folds for each entry in genomic coordinates
for entry in Coord:
	entryinfo = entry.rstrip().rsplit(sep = '\t')
	start = int(entryinfo[4])
	end = int(entryinfo[5])

	#Initialize parameters
	E40SummC = 0
	E55SummC = 0
	E40NmC = 0
	E55NmC = 0
	E40ICM.seek(0, 0)
	E55epi.seek(0, 0)
	E40ICMmC = E40ICM.readline().rsplit(sep = '\t')
	E55epimC = E55epi.readline().rsplit(sep = '\t')

	#Confirm the starting methylation site
	while E40ICMmC[0] != 'chr6' or int(E40ICMmC[1]) < start:
		E40ICMmC = E40ICM.readline().rsplit(sep = '\t')
	
	#Calculate mean methylation level of all sites within the range of certain gene
	while int(E40ICMmC[1]) <= end:
		E40SummC += float(E40ICMmC[3])
		E40NmC += 1
		E40ICMmC = E40ICM.readline().rsplit(sep = '\t')
	try:
		E40mean = E40SummC/E40NmC
	except:
		E40mean = 0
	while E55epimC[0] != 'chr6' or int(E55epimC[1]) < start:
		E55epimC = E55epi.readline().rsplit(sep = '\t')
	while int(E55epimC[1]) <= end:
		E55SummC += float(E55epimC[3])
		E55NmC += 1
		E55epimC = E55epi.readline().rsplit(sep = '\t')
	try:
		E55mean = E55SummC/E55NmC
	except:
		E55mean = 0
	if E40mean == 0:
		ratio = 'NA'
	else:
		ratio = E55mean/E40mean

	#Write the fold change into output file
	result = str(ratio) + '\t' + entry
	output.write(result)

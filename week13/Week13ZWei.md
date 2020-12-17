# Assignment 10: Metagenomics

**Zelin Wei**

17 Dec., 2020

## Download and inspect data

Download the data:

	mkdir Week13
	cd Week13
	wget https://bx.bio.jhu.edu/data/cmdb-lab/week13_data.tar
	tar -xvf week13_data.tar
	
Inspect `week13_data/KRAKEN/assembly.kraken`:

	less -S week13_data/KRAKEN/assembly.kraken
	
	>NODE_2_length_556123_cov_2361.439230    root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecalis;Enterococcus faecalis V583
	>NODE_11_length_278925_cov_118.155370    root;cellular organisms;Bacteria;Terrabacteria group;Actinobacteria;Actinobacteria;Propionibacteriales;Propionibacteriaceae;Cutibacterium;Cutibacterium avidum;Cutibacterium avidum 44067
	...
	
We can see that this file gives the taxonomy of each node in `assembly.fasta`.

## Visualize taxonomic profile

I plan to use `KronaTools` to visualize the taxonomic profile. However, `KronaTools` require that the input file is a tab-separated chart with first column being the total number of reads from one specie and the rest of the columns being the taxonomic information. To transform `.kraken` files to `KronaTools` acceptable form, I write a `python` script [`KrakenToKrona.py`](KrakenToKrona.py), and run:

	python KrakenToKrona.py week13_data/KRAKEN/SRR*
	ktImportText -o KronaCharts.html week13_data/KRAKEN/*.kronachart
	
Then I get the [pie charts](KronaCharts.html) as below:

<center><embed src=KronaCharts.html width=800 height=800></center>

--

### Question 1: Describe the trends in gut microbiota.

It can be seen that three main clades of micrbiomes exist at birth: the biggest clade is ***Lactobacillales*** being the overwhelmingly largest group making up 63% of total microbiomes (*E.faecalis* is the major specie for this clade); two other major clades are ***Bacillales*** and ***Proprionibacteriales***, with *S. epidermis* and *C. avidum* as their major species, respectively.

The percentage of ***Lactobacillales*** grows to 92% during the first day after birth, compressing ***Bacillales*** to only 7%. ***Propionibacteriales*** is further compressed to less than 1%.

From day 2 to day 5 after birth, the percentage of ***Bacillales*** gradually increases to 19%. On day 5 the percentage of ***Propionibacteriales*** dramatically rises to 3%.

Through day 5 to day 7, the percentage of ***Bacillales*** and ***Propionibacteriales*** continues to increase, reaching 10% and 28% respectively in day 7 after birth.

In general, all major clades and species exist shortly after birth. ***Lactobacillales*** overwhelms other clades during early stages, while ***Bacillales*** and ***Propionibacteriales*** gradually increase in their percentage. The diversity of gut microbiomes continuously increase in genral during the first week after birth, but the relative abundance of major species in each clades remains substantially unchanged.

--

## Metagenomic data binning

Then we need to group contigs from the same or closely related species together, i.e. binning the contigs.

--

### Question 2: What metrices can be used for binning?

1. Tetranucleotide frequency. The frequency of the four nucleotides should be similar across the genome of one specie, so contigs with similar tetranucleotide frequencies tend to be included in the same group. Similarly, GC content can also be used.
2. Co-abundance pattern. As described by [Sharon, *et al*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530670/), contigs sharing similar trends in abundance through time or cross samples should fall in the same group.
3. Coverage and coverage variance of the contigs.
4. (My guess) Preference in codon usage might be another parameter for binning. Closely related contigs theoretically should show similar preference in codon usage.

Combining different metrices would give more accurate binning results.

--
First I need to align reads to `assembly.fastq`. I use `bwa` to complete the alignment:

	bwa index -p assembly -a is week13_data/assembly.fasta 
	#Since assembly.fastq is not too big (< 2 Gbp), I use IS algorithm here.
	
	#Map paired reads to assembly.fastq.
	for file in week13_data/READS/*_1.fastq
	do
		temp=${file##*/}
		bwa mem -t 4 assembly $file ${file%%_1*}_2.fastq > ${temp%%_*}.mem
	done
	
	#Generate sorted .bam files as required by metabat2.
	for mem in *.mem
	do
		samtools sort -o ${mem%%.*}.bam -O BAM $mem
	done

Then, use `metabat2` to process sorted `.bam` files:

	runMetaBat.sh week13_data/assembly.fasta *.bam

--

### Question 3: 
#### A) How many bins?

I get 6 bins.

#### B) What percentage of contigs does each bin represent?

Use the following code to see how many contigs are included in each bin:

	cd assembly.fasta.metabat-bins-20201218_012439
	for bin in bin*
	do
		echo $bin
		grep -o '>NODE' $bin | wc -l
	done
	
	#Results:
	>bin.1.fa
      55
	>bin.2.fa
      78
	>bin.3.fa
       8
	>bin.4.fa
      37
	>bin.5.fa
      13
	>bin.6.fa
       6

There are in total 4103 contigs in `assembly.fasta`, so only 4.8% of the contigs are included by 6 bins. However, if I use `ls -l` to get the size of `bin.fa` files and `assembly.fasta` files, I find that nearly 98% of the sequences are covered by these bins:

	ls -l assembly.fasta.metabat-bins-20201218_012439
	
	>total 26224
		-rw-r--r--@ 1 zelinwei  staff  2752195 12 18 01:24 bin.1.fa
		-rw-r--r--  1 zelinwei  staff  2292366 12 18 01:24 bin.2.fa
		-rw-r--r--@ 1 zelinwei  staff  1683938 12 18 01:24 bin.3.fa
		-rw-r--r--  1 zelinwei  staff  1249747 12 18 01:24 bin.4.fa
		-rw-r--r--  1 zelinwei  staff  2525551 12 18 01:24 bin.5.fa
		-rw-r--r--@ 1 zelinwei  staff  2910800 12 18 01:24 bin.6.fa
	
	ls -l week13_data/assembly.fasta
	
	>	-rwxr-xr-x@ 1 zelinwei  staff  38856945 12  2  2017 week13_data/assembly.fasta

This indicates that these 6 bins include most of the major large contigs while small shattered contigs remain un-grouped.

#### C) Do bin sizes look right?

Prokaryotic genomes are about the size of several Mbp (*K. aerogenes* is 5.19 Mbp, *E. coli* is 4.6 Mbp). These bin sizes are of the same magnitude, so they seem reliable.

#### D) How to evaluate bin quality?

We can directly compare bin sequences with genomic sequences in databases. Good match would suggest high binning quality.

--

## Taxonomy of putative genomes

Download binning results:

	wget https://bx.bio.jhu.edu/data/cmdb-lab/bins.tar
	tar -xvf bins.tar

Map nodes in each bins to `assembly.kraken`:

	for bin in bins/*
	do
		grep NODE $bin > temp
		touch ${bin%.*}.kraken
		while read line
		do
			grep ${line#*>} week13_data/KRAKEN/assembly.kraken >> ${bin%.*}.kraken
		done < temp
		rm temp
	done

Use [`KrakenToKrona.py`](KrakenToKrona.py) to parse the output files and visualize by `KronaTools`:

	python KrakenToKrona.py bins/*.kraken
	ktImportText -o BinOrgs.html bins/*.kronachart

--

### Question 4:

#### A) Organism for each bin?

As we can see in [`BinOrgs.html`](BinOrgs.html) below:

<center><embed src=BinOrgs.html width=800 height=800>

|Bin Number|Representative specie|
|:----:|:----:|
|1|*Staphylococcus haemolyticus*|
|2|*Leuconostoc citreum*|
|3|*Staphylococcus lugdunensis*|
|4|*Enterococcus faecalis*|
|5|*Cutibacterium avidum*|
|6|*Staphylococcus epidermidis*|
|7|*Staphylococcus aureus subsp*|
|8|*Anaerococcus prevotii*|

</center>

#### B) More robust taxonomy method?

Use `BLAST` to align bin sequences to bacteria genomes in databases. `BLAST` will give scores so this approach might be more quantitative.

## Heatmap of bin abundance across time

Download bin abundance data:

	wget https://bx.bio.jhu.edu/data/cmdb-lab/abundance_table.tab

I use a `Python` script [`BinAbundanceHeatmap.py`](BinAbundanceHeatmap.py) to draw the heatmap:

	python BinAbundanceHeatmap.py
	
Please note that the data are log-transformed due to the huge variance among bins. The [heatmap](BinAbundance.png) is shown below:

<center><img src=BinAbundance.png></center>

--

### Question 5: Compare the heatmap with `Krona` pie charts.

The three main clades shown in `Krona` charts (represented by *E.faecalis*, *S. epidermis*, and *C. avidum*, respectively) also have high abundance in this heatmap. The heatmap also reflects the sharp increase of *C. avidum* percentage around 'SRR492193', and the continuous increase of *S. epidermis* percentage.

However, heatmap also includes more earlier samples ('SRR492065'-'SRR492182'), and we can see that *S. epidermis* and *C. avidum* abundance actually has a high plateau period during early stages--this is not presented by previous charts.

In addition, heatmap shows more details. For example, *L. citreum* abundance seems to reach a small plateau between 'SRR492184' and 'SRR492187', and then decreases to the bottom around 'SRR492193'. Such pattern seems to negatively correlates with *C. avidum*. Correlations between the abundance of different species might reflect inter-specie competition or mutualism relationship.
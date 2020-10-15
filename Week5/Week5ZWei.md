# Week 5 Assignments: ChIP-seq and Motif Finding

**Zelin Wei**

12 Oct., 2020

Fisrt, download the data needed:

	wget http://67.207.142.119/outgoing/g1e.tar.xz
	tar -xf g1e.tar.xz

## Part 1: ChIP-seq

### Mapping reads

Download reference file:

	wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr19.fa.gz
	gunzip chr19.fa.gz

Generate index file using:

	bowtie2-build chr19.fa chr19idx
	mkdir chr19idx
	mv chr19idx* chr19idx/

Then map each sample to the reference genome:

	for file in *.fastq
	do
		bowtie2 -x chr19idx/chr19idx -U $file -S ${file%.*}.sam -p 6
	done

### Calling peaks

Call peaks for each condition and write the output files into different directories:

	mkdir G1E
	mkdir ER4
	macs2 callpeak -t CTCF_G1E.sam -c input_G1E.sam -f SAM --outdir G1E -n G1E -s 36
	macs2 callpeak -t CTCF_ER4.sam -c input_ER4.sam -f SAM --outdir ER4 -n ER4 -s 36
	
The `narrowPeak` files are: [`G1E_peaks.narrowPeak`](G1E/G1E_peaks.narrowPeak) and [`ER4_peaks.narrowPeak`](ER4/ER4_peaks.narrowPeak).
	
### Differential binding

Use `bedtools` to spot peaks with differential binding under different conditions:

	bedtools intersect -a G1E/G1E_peaks.narrowPeak -b ER4/ER4_peaks.narrowPeak -v > Losspeaks.bed
	# Find lost peaks.
	bedtools intersect -a ER4/ER4_peaks.narrowPeak -b G1E/G1E_peaks.narrowPeak -v > Gainpeaks.bed
	# Find gained peaks.

The output files are [`Losspeaks.bed`](Losspeaks.bed) and [`Gainpeaks.bed`](Gainpeaks.bed).
	
### Feature overlapping

Download annotation file:

	wget https://raw.githubusercontent.com/bxlab/qbb2020/master/week5/Mus_musculus.GRCm38.94_features.bed
	
Overlap the differential binding peaks and the annotation file:

	bedtools intersect -a Losspeaks.bed -b Mus_musculus.GRCm38.94_features.bed -wa -wb > LosspeaksAnnot.bed
	bedtools intersect -a Gainpeaks.bed -b Mus_musculus.GRCm38.94_features.bed -wa -wb > GainpeaksAnnot.bed
	
Count the number of each feature:

	cut -f 14 LosspeaksAnnot.bed | sort | uniq -c
	cut -f 14 GainpeaksAnnot.bed | sort | uniq -c
	
The results are:

---

<center>
**Feature counts for peaks lost in ER4-induced differentiation**

|Feature|Peak count number|
|:----:|:----:|
|Exon|10|
|Intron|30|
|Promoter|2|

--

**Feature counts for peaks gained in ER4-induced differentiation**

|Feature|Peak count number|
|:----:|:----:|
|Exon|31|
|Intron|72|
|Promoter|23|

</center>

---

### Plotting

Count features for each condition (G1E or ER4):

	for sample in G1E ER4
	do
		bedtools intersect -a ${sample}/${sample}_peaks.narrowPeak -b Mus_musculus.GRCm38.94_features.bed -wa -wb | cut -f 14 | sort | uniq -c > ${sample}FeatureCount.csv
	done
	
Count differential binding peak numbers:

	wc -l Losspeaks.bed Gainpeaks.bed > Diffpeaks.csv
	
Write a `Python` script [`PeakPlot.py`](PeakPlot.py) to draw the plots:

	python PeakPlot.py
	
The plots are shown below:

<center>

<img src=PeakPlot.png>

</center>

## Part 2: Motif discovery

Download motif databases:

	wget http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.19.tgz
	mkdir MotifDB
	tar xvf motif_databases.12.19.tgz.html -C MotifDB/
	
Prepare sequence file for the top 100 CTCF peaks in ER4-treated cells:

	sort -k7,7nr ER4/ER4_peaks.narrowPeak | head -100 > ER4/ER4TopPeaks.bed
	bedtools getfasta -fi chr19.fa -bed ER4/ER4TopPeaks.bed -fo ER4TopPeaksSeq.fa
	
Use `meme-chip` to search for potential motifs:

	meme-chip -meme-maxw 20 ER4TopPeaksSeq.fa
	
Compare the predicted motifs with databases:

	tomtom memechip_out/meme_out/meme.txt MotifDB/motif_databases/JASPAR/JASPAR_CORE_2016.meme
	
From `tomtom.tsv`, motif 'RCCACYAGGKGGCRCYVKDG' has the lowest p-value:

	sort -k4,4g tomtom_out/tomtom.tsv | head

Therefore I drew a sequence logo for this motif:

	ceqlogo -iRCCACYAGGKGGCRCYVKDG memechip_out/meme_out/meme.txt -o TopMotifLogo.eps -f EPS
	epstopdf TopMotifLogo.eps
	
Shown below is the sequence logo for 'RCCACYAGGKGGCRCYVKDG' (aligned to 'MA0139.1' in `JASPAR_CORE_2016.meme`):

<center>

<img src=TopMotifLogo.pdf>

</center>
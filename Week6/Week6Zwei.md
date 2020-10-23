# Week 6 Assignment: DNA Methylation

**Zelin Wei**

20 Oct., 2020

First build up the environment:

	conda create -n week6 -c bioconda -c anaconda fastqc bismark samtools bowtie2 igv
	conda activate week6
	
## Download data

Download the data. Based on the [GEO Accession viewer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76505), <font color=#FF0000>`SRR3083926`</font> is **STEM-seq E4.0ICM rep1**, and <font color=#FF0000>`SRR3083929`</font> is **STEM-seq E5.5Epi rep1**. They are both paired reads.

Run `fastqc` on `SRR3083926_1.chr6.fastq`, and the report shows that the C percentage is very low in the `Per base sequence content` category. Since this is a bisulfite sequencing data, this is likely the result of 5mC-to-U transition.

## Bisulfite mapping

### Prepare reference genome

Download chromosome 6 in mm10:

	wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr6.fa.gz
	gunzip chr6.fa.gz
	mkdir mm10chr6
	mv chr6.fa mm10chr6/
	
Create index files for the reference genome:

	bismark_genome_preparation mm10chr6/

### Map bisulfite reads to genome

Map the reads using `bismark`:

	bismark mm10chr6 -1 methylation_fastq/SRR3083926_1.chr6.fastq,methylation_fastq/SRR3083929_1.chr6.fastq -2 methylation_fastq/SRR3083926_2.chr6.fastq,methylation_fastq/SRR3083929_2.chr6.fastq
	
To visualize the data in `IGV`, sort the `.bam` files:

	samtools sort SRR3083926_1.chr6_bismark_bt2_pe.bam > SRR3083926_1.chr6_bismark_bt2_pe_sorted.bam
	samtools sort SRR3083929_1.chr6_bismark_bt2_pe.bam > SRR3083929_1.chr6_bismark_bt2_pe_sorted.bam
	samtools index SRR3083926_1.chr6_bismark_bt2_pe_sorted.bam
	samtools index SRR3083929_1.chr6_bismark_bt2_pe_sorted.bam
	
### Extract methylated bases

Extract methylated C in three contexts: CpG, CHG, and CHH, and merge these data into a `.bedgraph` file:

	bismark_methylation_extractor --bedGraph --comprehensive *pe.bam
	gunzip *.gz
	
View the `.bedgraph` file in `IGV`. Below is a snapshot for **chr6: 50,000,000-60,000,000**:

<center><img src=ExtractMethylation.png></center>

### Calculate methylation fold change

Download gene coordinates in **chr6: 50,000,000-60,000,000**:

	wget https://drive.google.com/file/d/1vbpa5EehoG3Yt91lNJqiiSsvtkkA3d0x/view?usp=sharing
	
Write a `Python` script [`MethylRatio.py`](MethylRatio.py) to calculate average methylation level fold change on each gene and then sort the results. Fold change of each genomic coordinate entry is added to the first column of the original coordinate `.bed` file:

	python MethylRatio.py
	sort -k1,1nr MethylRatio.bed > MethylRatioSorted.bed
	
In [`MethylRatioSorted.bed`](MethylRatioSorted.bed) we can see that the three entries with a methylation fold change larger than 15 are **NR_046179**, **NM_011153**, and **NR_110440**. Genes included in these entries are: ***6430584L05Rik***, ***Ppp1r17***, ***cmpl***, and ***Halr1***.
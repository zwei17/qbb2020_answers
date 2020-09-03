#Prepare Conda environment
conda create --name day3
conda activate day3
conda install samtools=1.1
which samtools

#Prepare workspace
cd
cd qbb2020_answers
mkdir day3_hw
cp ../qbb2020/data/K*.bed day3_hw 
cp ../qbb2020/data/BDGP6.Ensembl.81.gtf day3_hw 
cd day3_hw
wget https://bx.bio.jhu.edu/track-hubs/cmdb/SRR072893.bam
samtools index SRR072893.bam
ls -l > directory.txt

#Visualize regions flanking sxl
conda install igv
igv

#Count features in BDGP.Ensembl.81.gtf
#cut -f 3 BDGP6.Ensembl.81.gtf | sort | uniq -c | tail -8 > features.txt
cut -f 3 BDGP6.Ensembl.81.gtf | sort | uniq -c | grep -v '#' > features.txt #Revised based on the presentation.


#Count total interval numbers
for fname in *.bed
do
	#echo Total interval number of $filename  >> ${filename/.bed/.info}
	#wc -l $filename >> ${filename/.bed/.info}
	cut -f 1 fname | sort | uniq -c > ${filename/.bed/.info} #Revised based on the presentation.
done

#Report first 10 intervals on 2L
grep 2L K4me3.bed | sort -k 2 -n | head -10 > 2L_first10intervals.txt

#Install bedtools for advanced questions
conda install bedtools

#Sort bed files
for fname in *.bed
do
	sort -k1,1 -k2,2n $filename > ${filename/./.sorted.}
done

#Calculate Jaccard correlation among bed files
echo '*** K4 vs K9' >> jaccard.txt
bedtools jaccard -a K4me3.sorted.bed -b K9me3.sorted.bed >> jaccard.txt
echo '*** K4 vs K27' >> jaccard.txt
bedtools jaccard -a K4me3.sorted.bed -b K27me3.sorted.bed >> jaccard.txt
echo '*** K9 vs K27' >> jaccard.txt
bedtools jaccard -a K9me3.sorted.bed -b K27me3.sorted.bed >> jaccard.txt

#Annotate ChIP seq data
cp ../../qbb2020/data/genes.bed .
bedtools closest -a K4me3.sorted.bed -b genes.bed -d -t first | cut -f -1,2,3,7,8 > K4me3_annotated.bed
tail -10 K4me3_annotated.bed > K4me3_annotated_last10.bed

#Find intervals with highest expression
bedtools intersect -a K4me3_annotated.bed -b SRR072893.bam -c > K4me3_counts.bed
sort -k6 -n K4me3_counts.bed | tail -10 | cut -f -1,2,3,4,6 > K4me3_highexp.txt
#Check act5c and rps3a in igv, images committed.
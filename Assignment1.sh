#usr/bin/awk
##Copy a directory providing the paired-end sequence data and the sample details.
cd
cp -r /localdisk/data/BPSM/Assignment1/fastq rawdata
#Check files containing raw sequence data and sample information.
cd ./rawdata
ls -lh ./fastq
cat ./fastq/fqfiles
#perform a quality check on the paired-end raw sequence data.
mkdir result_1
gzip -d ./fastq/2*fq.gz
fastqc -f fastq -o result_1 ./fastq/2*fq
#align the read pairs to the Trypanosoma brucei brucei genome after assessing the number and #quality of the raw sequence data based on the output of the last step.
cp -r /localdisk/data/BPSM/Assignment1/Tbb_genome .
ls ./Tbb_genome
gzip -d ./Tbb_genome/Tb927_genome.fasta.gz
#build a bowtie index for alignment.
bowtie2-build ./Tbb_genome/Tb927_genome.fasta Tbbgenome
#align paired-end example reads
touch alignment1.sam
touch alignment2.sam
touch alignment3.sam
touch alignment4.sam
touch alignment5.sam
touch alignment6.sam
bowtie2 --fast -x Tbbgenome -1 ./fastq/216_L8_1.fq -2 ./fastq/216_L8_2.fq -S alignment1.sam
bowtie2 --fast -x Tbbgenome -1 ./fastq/218_L8_1.fq -2 ./fastq/218_L8_2.fq -S alignment2.sam
bowtie2 --fast -x Tbbgenome -1 ./fastq/219_L8_1.fq -2 ./fastq/219_L8_2.fq -S alignment3.sam
bowtie2 --fast -x Tbbgenome -1 ./fastq/220_L8_1.fq -2 ./fastq/220_L8_2.fq -S alignment4.sam
bowtie2 --fast -x Tbbgenome -1 ./fastq/221_L8_1.fq -2 ./fastq/221_L8_2.fq -S alignment5.sam
bowtie2 --fast -x Tbbgenome -1 ./fastq/222_L8_1.fq -2 ./fastq/222_L8_2.fq -S alignment6.sam
#convert the output to indexed "bam" format.
samtools view -bS alignment1.sam > alignment1.bam | samtools sort -o alignment1.sorted.bam
samtools view -bS alignment2.sam > alignment2.bam | samtools sort -o alignment2.sorted.bam
samtools view -bS alignment3.sam > alignment3.bam | samtools sort -o alignment3.sorted.bam
samtools view -bS alignment4.sam > alignment4.bam | samtools sort -o alignment4.sorted.bam
samtools view -bS alignment5.sam > alignment5.bam | samtools sort -o alignment5.sorted.bam
samtools view -bS alignment6.sam > alignment6.bam | samtools sort -o alignment6.sorted.bam
samtools index alignment1.sorted.bam
samtools index alignment2.sorted.bam
samtools index alignment3.sorted.bam
samtools index alignment4.sorted.bam
samtools index alignment5.sorted.bam
samtools index alignment6.sorted.bam
#count the number of reads that align to the regions of the genome that code for genes.
cp /localdisk/data/BPSM/Assignment1/Tbbgenes.bed
mkdir result_2
bedtools multicov -bams alignment1.sorted.bam alignment2.sorted.bam alignment3.sorted.bam -bed ./Tbbgenes.bed > ./result_2/slender
bedtools multicov -bams alignment4.sorted.bam alignment5.sorted.bam alignment6.sorted.bam -bed ./Tbbgenes.bed > ./result_2/stumpy
#count the statistical mean of per gene for each group.
head -n6 ./result_2/slender
cut -f 4,7,8,9 ./result_2/slender > ./slender_samples.txt
cut -f 4,7,8,9 ./result_2/stumpy > ./stumpy_samples.txt
awk 'NF==1 {next} {T = 0; for(i=2; i <= NF; i++) T += $i; T /= (NF-1); print $1, T}' ./slender_samples.txt
awk 'NF==1 {next} {T = 0; for(i=2; i <= NF; i++) T += $i; T /= (NF-1); print $1, T}' ./stumpy_samples.txt
paste ./slender_samples.txt ./stumpy_samples.txt > mean.txt
#the left means of per gene are for slender group
#the right means of per gene are for stumpy group




 
 
 





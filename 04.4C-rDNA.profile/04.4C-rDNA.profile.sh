#!/bin/sh
# script to get and process HEK293T 4C-rDNA reads
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
#
# Output: 1. 4C.HEK293.bedGraph, 4C.HEK293.bw    genome-wide hg38 HEK293 4C-rDNA profile for genome browsers and/or profile charting
# Dependency tools:
# 1. NCBI SRA Toolkit     https://github.com/ncbi/sra-tools
# 2. cutadapt             https://cutadapt.readthedocs.io
# 4. bwa                  https://github.com/lh3/bwa
# 5. samtools             http://www.htslib.org
# 6. bedtools 2.29.1      https://bedtools.readthedocs.io/en/latest/
# 7. UCSC Kent's tools    https://hgdownload.cse.ucsc.edu/admin/exe/
# 8. DFAM database hg38   https://www.dfam.org/releases/Dfam_3.4/annotations/hg38/hg38_dfam.nrph.hits.gz
#
# Get data from SRA Archive
# SRX4900048:GSM3434713 4C-rDNA HEK293T replica 1
fastq-dump --split-files --gzip SRR8072070
# SRX4900049:GSM3434714 4C-rDNA HEK293T replica 1
fastq-dump --split-files --gzip SRR8072071
mv SRR8072070_1.fastq.gz rep1.fastq.gz
mv SRR8072071_1.fastq.gz rep2.fastq.gz
# Quality and adapters trimming for replica 1
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=4C_notail.fastq.gz \
-g GCCTAAGCCTGCTGAGAACTTTC -g CAGCATTCTGTAGGGAGATCAAATC -a GAAAGTTCTCAGCAGGCTTAGGC -a GATTTGATCTCCCTACAGAATGCTG \
-o 4C_tail.fastq.gz rep1.fastq.gz
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=cutdebug.fastq.gz \
-g TCTTTGAAAAAAATCCCAGAAGTGGT -g AAGTCCAGAAATCAACTCGCCAGT -a ACTGGCGAGTTGATTTCTGGACTT -a ACCACTTCTGGGATTTTTTTCAAAGA \
-o 4C_head.fastq.gz 4C_notail.fastq.gz
cat 4C_tail.fastq.gz 4C_head.fastq.gz > 4C_trimmed.rep1.fastq.gz
rm -f 4C_tail.fastq.gz 4C_head.fastq.gz 4C_notail.fastq.gz cutdebug.fastq.gz
# Quality and adapters trimming for replica 2
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=4C_notail.fastq.gz \
-g GCCTAAGCCTGCTGAGAACTTTC -g CAGCATTCTGTAGGGAGATCAAATC -a GAAAGTTCTCAGCAGGCTTAGGC -a GATTTGATCTCCCTACAGAATGCTG \
-o 4C_tail.fastq.gz rep2.fastq.gz
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=cutdebug.fastq.gz \
-g TCTTTGAAAAAAATCCCAGAAGTGGT -g AAGTCCAGAAATCAACTCGCCAGT -a ACTGGCGAGTTGATTTCTGGACTT -a ACCACTTCTGGGATTTTTTTCAAAGA \
-o 4C_head.fastq.gz 4C_notail.fastq.gz
cat 4C_tail.fastq.gz 4C_head.fastq.gz > 4C_trimmed.rep2.fastq.gz
rm -f 4C_tail.fastq.gz 4C_head.fastq.gz 4C_notail.fastq.gz cutdebug.fastq.gz
# alignment to hg38 replica 1
bwa mem -t 22 /usr/local/genomes/hg38.mfa 4C_trimmed.rep1.fastq.gz | samtools view -bS -F 4 -o hits.bam
samtools sort -@ 10 hits30-200.bam -O BAM -o hits.rep1.bam
# alignment to hg38 replica 2
bwa mem -t 22 /usr/local/genomes/hg38.mfa 4C_trimmed.rep2.fastq.gz | samtools view -bS -F 4 -o hits.bam
samtools sort -@ 10 hits30-200.bam -O BAM -o hits.rep2.bam
rm -f hits.bam
# Convert BAM files to bedGraph
genomeCoverageBed -bg -ibam hits.rep1.bam > rep1.bedGraph
genomeCoverageBed -bg -ibam hits.rep2.bam > rep2.bedGraph
# get DFAM database and convert it to bed format
wget https://www.dfam.org/releases/Dfam_3.4/annotations/hg38/hg38_dfam.nrph.hits.gz
./lib/dfam2bed hg38_dfam.nrph.hits.gz
mv hg38_dfam.nrph.hits.bed hg38_dfam.bed
# Building intersection and DFAM removing
sort -k1,1 -k2,2n -k3,3n rep1.bedGraph -o rep1.bedGraph
sort -k1,1 -k2,2n -k3,3n rep2.bedGraph -o rep2.bedGraph
bedtools intersect -sorted -v -f 1.0 -a rep1.bedGraph -b hg38_dfam.bed  > 4C_hg38_nodfam.rep1.bedGraph
bedtools intersect -sorted -v -f 1.0 -a rep2.bedGraph -b hg38_dfam.bed  > 4C_hg38_nodfam.rep2.bedGraph
bedtools intersect -a 4C_hg38_nodfam.rep1.bedGraph -b 4C_hg38_nodfam.rep2.bedGraph -wb > intersect_reps.txt
./lib/overlist2bed_mean.pl intersect_reps.txt
sort -k1,1 -k2,2n -k3,3n intersect_reps.bed -o 4C.HEK293.bedGraph
rm -f intersect* hits*.bam
bedGraphToBigWig 4C.HEK293.bedGraph ~/bin/hg38.chr.sizes 4C.HEK293.bw

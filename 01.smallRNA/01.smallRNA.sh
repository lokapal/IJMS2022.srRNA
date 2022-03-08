#!/bin/sh
# script to get and process HEK293T smallRNA reads
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Output: 1. table_srRNA.txt the file with the genome-wide hg38 srRNA mappings: chromosome coordinates, alignment length, coverage, reads, sequence
#         2. srRNA_density.bedGraph, srRNA_density.bw            genome-wide hg38 srRNA profile for genome browsers and/or profile charting
#            srRNA_density_frw.bedGraph, srRNA_density_frw.bw    genome-wide hg38 forward strand srRNA profile for genome browsers and/or profile charting
#            srRNA_density_rev.bedGraph, srRNA_density_rev.bw    genome-wide hg38 reverse strand srRNA profile for genome browsers and/or profile charting
#
# Dependency tools:
# 1. NCBI SRA Toolkit     https://github.com/ncbi/sra-tools
# 2. Trimmomatic 0.36     http://www.usadellab.org/cms/?page=trimmomatic
# 3. BBTools              https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
# 4. bowtie2 2.3.4.1      http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# 5. samtools             http://www.htslib.org
# 6. bedtools 2.29.1      https://bedtools.readthedocs.io/en/latest/
# 7. UCSC Kent's tools    https://hgdownload.cse.ucsc.edu/admin/exe/
# 8. Perl with BIO::DB::Fasta BioPerl library
#
#SRX157657 GSM955512 smallRNA RNA-Seq
# Get data from SRA Archive
fastq-dump --split-files --gzip SRR518497
# Quality trimming
trimmomatic SE -threads 10 SRR518497_1.fastq.gz smallRNA.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
# Exact copies deduplication
dedupe.sh in=smallRNA.filtered.fastq.gz out=unique_smallRNA.fastq.gz ac=f qtrim=rl usejni=t overwrite=t
#Alignment of non-unique sequences to H.sapiens rDNA Genbank U13369 to detect srRNA amongst all small RNA
bowtie2 --end-to-end --very-sensitive -p 22 -x /usr/local/genomes/hs_rDNA -U smallRNA.filtered.fastq.gz --no-unal | samtools view -bS -F 4 -o hits.bam
samtools sort -@ 20 hits.bam -O BAM -o hits_sorted.bam
mv hits_sorted.bam hits.bam
#Conversion aligned reads back to FastQ format
bamToFastq -i hits.bam -fq end2end.fastq
gzip end2end.fastq
#Alignment of srRNA to H.sapiens hg38 genome
bowtie2 --end-to-end --sensitive -p 22 -x /usr/local/genomes/hg38 -U end2end.fastq.gz --no-unal | samtools view -bS -F 4 -o hits.bam
samtools sort -@ 20 hits.bam -O BAM -o hits_sorted.bam
mv hits_sorted.bam hits.bam
samtools mpileup -f /usr/local/genomes/hg38.mfa hits.bam -o pileup.txt
samtools view -@ 20 -O SAM hits.bam -o hits.sam
samtools index -@ 20 hits.bam
sort -k1,1 -k2,2n pileup.txt | bgzip -c > compressed.pileup.gz
tabix -s 1 -b 2 -e 2 compressed.pileup.gz
perl ./lib/makeTablePileup.pl pileup.txt > table.txt
perl ./lib/addNumOfReads.pl hits.sam table.txt > tableA.txt
perl ./lib/addSubseq.pl /usr/local/genomes/hg38.mfa tableA.txt > table_srRNA.txt
rm -f compressed.* hits.sam table.txt tableA.txt pileup.*
#Genome-wide srRNA profiles creation for the forward and reverse strands
# Dividing forward and reverse strands
samtools view -F 0x10 hits.bam -o forward.bam
samtools view -f 0x10 hits.bam -o reverse.bam
#convert to bedGraph
genomeCoverageBed -bg -ibam hits.bam > srRNA_density.bedGraph
genomeCoverageBed -bg -ibam forward.bam > srRNA_density_frw.bedGraph
genomeCoverageBed -bg -ibam reverse.bam > srRNA_density_rev.bedGraph
#sort
bedSort srRNA_density.bedGraph srRNA_density.bedGraph
bedSort srRNA_density_frw.bedGraph srRNA_density_frw.bedGraph
bedSort srRNA_density_rev.bedGraph srRNA_density_rev.bedGraph 
#convert
bedGraphToBigWig srRNA_density.bedGraph ~/bin/hg38.chr.sizes srRNA_density.bw
bedGraphToBigWig srRNA_density_frw.bedGraph ~/bin/hg38.chr.sizes srRNA_density_frw.bw
bedGraphToBigWig srRNA_density_rev.bedGraph ~/bin/hg38.chr.sizes srRNA_density_rev.bw

#!/bin/bash
pwd
inputs=""
aligns=""
isinputs=0
if [ -d MACS2 ]
then
exit
fi
if [ -e *Rep1.fastq.gz ]
then
echo Processing Rep1...
pattern="*Rep1.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 22 $infile Rep1.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=Rep1.filtered.fastq.gz out=unique1.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 unique1.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam hits1.bam
aligns+=" hits1.bam"
fi
if [ -e *Rep2.fastq.gz ] 
then
echo Processing Rep2...
pattern="*Rep2.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 10 $infile Rep2.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=Rep2.filtered.fastq.gz out=unique2.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 unique2.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam hits2.bam
aligns+=" hits2.bam"
fi
if [ -e *Rep3.fastq.gz ] 
then
echo Processing Rep3...
pattern="*Rep3.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 22 $infile Rep3.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=Rep3.filtered.fastq.gz out=unique3.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 unique3.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam hits3.bam
aligns+=" hits3.bam"
fi
if [ -e *Rep4.fastq.gz ] 
then
echo Processing Rep4...
pattern="*Rep4.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 22 $infile Rep4.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=Rep4.filtered.fastq.gz out=unique4.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 unique4.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam hits4.bam
aligns+=" hits4.bam"
fi
if [ -e *input.fastq.gz ]
then
echo Processing input...
pattern="*input.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 22 $infile input.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=input.filtered.fastq.gz out=uniqueI.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 uniqueI.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam input.bam
#pigz -9 uniqueI.fastq
inputs+=" input.bam"
((isinputs++))
fi
if [ -e *input1.fastq.gz ]
then
echo Processing input1...
pattern="*input1.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 22 $infile input1.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=input1.filtered.fastq.gz out=uniqueI1.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 uniqueI1.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam input1.bam
inputs+=" input1.bam"
((isinputs++))
fi
if [ -e *input2.fastq.gz ]
then
echo Processing input2...
pattern="*input2.fastq.gz"
files=( $pattern )
infile="${files[0]}"
trimmomatic SE -threads 22 $infile input2.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
dedupe.sh in=input2.filtered.fastq.gz out=uniqueI2.fastq.gz ac=f qtrim=rl usejni=t overwrite=t # qin=33
bowtie --chunkmbs 1024 -q -p 22 -a --best -m 1 /usr/local/genomes/hg38 uniqueI2.fastq.gz -S | samtools view -bS -F 4 -o hits30-200.bam
samtools sort --threads 10 -o hits30-200sorted.bam hits30-200.bam
mv hits30-200sorted.bam input2.bam
inputs+=" input2.bam"
((isinputs++))
fi
# prepare commandstring for MACS2
currdir=${PWD##*/}
name="--name "
name+=$currdir
cstring="macs2 callpeak -f BAM --bdg --gsize hs --call-summits --outdir MACS2 "
cstring+=$name
cstring+=" -t "
cstring+=$aligns
if [ $isinputs -gt 0 ]
then
cstring+=" -c "
cstring+=$inputs
fi
echo $cstring
eval $cstring
rm -f *.bam

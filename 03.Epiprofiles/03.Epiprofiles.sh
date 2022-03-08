#!/bin/sh
# script to get and process ChIP-Seq 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# script to process all ChIP-Seq files like *Rep?.fastq.gz and *input?.fastq.gz
# Up to 4 IP data replicas and up to 2 inputs are supported, name wildcards should be:
# Inputs: IP data: *Rep1.fastq.gz
#                  *Rep2.fastq.gz
#                  *Rep3.fastq.gz
#                  *Rep4.fastq.gz
#     input data:  *input.fastq.gz
#                  *input1.fastq.gz
#                  *input2.fastq.gz
# 
# from raw ChIP-Seq reads to MACS2 https://github.com/taoliu/MACS peaks
# Dependency tools:
# 1. Trimmomatic 0.36   http://www.usadellab.org/cms/?page=trimmomatic
# 2. BBTools            https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
# 3. bowtie 1.3.1       http://bowtie-bio.sourceforge.net/index.shtml
# 4. samtools           http://www.htslib.org
# 5. MACS2              https://github.com/taoliu/MACS, https://doi.org/10.1186/gb-2008-9-9-r137
./1.chipseq.align.sh
./2.chipseq.diffprof.sh

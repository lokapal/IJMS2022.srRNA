#!/bin/bash
# script to create differential fold enrichment ChIP-Seq profiles/signals
# from MACS2 bedGraph intermediate files to BigWig files
# should be run from directory one level up from MACS2 processed files
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Dependency tools: MACS2              https://github.com/taoliu/MACS, https://doi.org/10.1186/gb-2008-9-9-r137
#                   UCSC Kent's tools  https://hgdownload.cse.ucsc.edu/admin/exe/
pwd
# prepare commandstring for MACS2
currdir=${PWD##*/}
cd MACS2
unpgz="gzip -d *.gz"
eval $unpgz
cstring="macs2 bdgcmp -t *treat_pileup.bdg  -c *control_lambda.bdg -o FE.bdg -m FE"
eval $cstring
pgz1="gzip *treat_pileup.bdg"
eval $pgz1
pgz2="gzip *control_lambda.bdg"
eval $pgz2
cstring3="sort -k1,1 -k2,2n FE.bdg > sorted.bdg"
eval $cstring3
rm FE.bdg
name=$currdir
name+=".bw"
cstring2="bedGraphToBigWig sorted.bdg ~/bin/hg38.chr.sizes "
cstring2+=$name
eval $cstring2
rm sorted.bdg
cd ..

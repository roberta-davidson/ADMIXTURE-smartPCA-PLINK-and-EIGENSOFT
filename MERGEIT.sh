#!/bin/bash

source activate eig

in1=$1
in2=$2
out=$3

mergeit -p <(echo "geno1: ${in1}.geno
snp1:  ${in1}.snp
ind1:  ${in1}.ind
geno2: ${in2}.geno
snp2:  ${in2}.snp
ind2:  ${in2}.ind
genooutfilename:   ${out}.geno
snpoutfilename:    ${out}.snp
indoutfilename:    ${out}.ind
outputformat:  EIGENSTRAT
docheck: YES
strandcheck:  YES
hashcheck: NO")

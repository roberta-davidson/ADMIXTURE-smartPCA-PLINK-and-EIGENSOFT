#!/bin/bash

source activate eig

in1=$1
out=$2

convertf -p <(echo "genotypename:	${in1}.geno
snpname:	${in1}.snp
indivname:	${in1}.ind
outputformat:	EIGENSTRAT 
genotypeoutname:	${out}.geno
snpoutname:	${out}.snp
indivoutname:	${out}.ind
poplistname:	poplist.txt") #poplist.txt is a list file of populations (column 3 of .ind file) to keep

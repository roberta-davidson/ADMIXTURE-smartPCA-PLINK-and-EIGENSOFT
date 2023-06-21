#!/bin/bash

source activate eig

# $1 = first file after "sbatch script.sh" when you submit on the command line

in1=$1
in2=$2
in3=$3

convertf -p <(echo "genotypename:	${in1}.geno
snpname:	${in1}.snp
indivname:	${in1}.ind
outputformat:	EIGENSTRAT
genotypeoutname:	${in2}.geno
snpoutname:	${in2}.snp
indivoutname:	${in2}.ind
#poplistname:	poplist.txt #poplist if you wish to subset by individuals at the same time
badsnpname:	${in3}") #list of SNPs to remove

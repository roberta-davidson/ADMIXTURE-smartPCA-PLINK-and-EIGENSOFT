#!/bin/bash

source activate eig

in1=$1

convertf -p <(echo "genotypename:	${in1}.geno
snpname:	${in1}.snp
indivname:	${in1}.ind
outputformat:	PACKEDPED
genotypeoutname:	${in1}.bed
snpoutname:	${in1}.bim
indivoutname:	${in1}.fam")

#!/bin/bash

source activate eig

in1=$1
in2=poplist.txt

smartpca -p <(echo "genotypename:	${in1}.geno
snpname:	${in1}.snp
indivname:	${in1}.ind
evecoutname:	${in1}.pca.evec.txt
evaloutname:	${in1}.pca.eval.txt
poplistname:	${in2}
missingmode:	NO
lsqproject:	YES
outliermode:	1
outlieroutname:	${in1}.pca.outlier.txt
#popsizelimit:	50
snpweightoutname:	${in1}.pca.snpweight.txt
numthreads:	4
newshrink:	YES
maxpops:	500
numoutevec:	10")

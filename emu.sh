#!/bin/bash

#requirements Python=3

conda activate emu

DATA=$1

emu --plink ${DATA} \
	--n_eig 10 \
	--threads 8 \
	--maf 0.01 \
	--out ${DATA} \
	--maf_save \
	--sites_save

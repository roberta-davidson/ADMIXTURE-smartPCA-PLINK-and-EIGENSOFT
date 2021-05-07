# ADMIXTURE-PCA-PLINK-EIGENSTRAT
Scripts and notes on how to work with PLINK files, EIGENSTRAT files, converting between them. |
How to run and plot ADMIXTURE and smartPCA. \

## Using PLINK
PLINK is in general very annoying, reccommend to manipulate data in VCF or in EIGENSTRAT formats where possible. \
There are many functions PLINK will do to your data by default, so find the flags necessary to turn off these functions. \
Some useful ones I use: \
`--keep-allele-order`	Use this EVERY SINGLE TIME you call a plink command, otherwise the order of Allele1 and Allele2 may (or probably will) be flipped in your data. \
`--allow-no-sex` 	PLINK will default to removing individuals that have unassigned sex, use this to force it to keep them. \
`--snps-only` 		Removes indels from your variant data and keeps only snps \
`--biallelic-only` Removes sites with 2+ alleles \
`--indiv-sort 0` PLINK default re-orders your data by individual name, this keeps them the same order as the `*.fam` file \
`--geno 0.9999`	Removes sites with greater that 0.9999 missing data, so effectively removes a locus with no data \
`--extract`/`--exclude` Extracts or exlcludes variants based on a .txt file list of all variant IDs

Keep the PLINK manual handy: https://www.cog-genomics.org/plink/1.9/

## Convert VCF to PLINK format
```
ml plink/1.90beta-4.4-21-May

plink \
  --vcf <in>.vcf.gz \
  --allow-no-sex \
  --keep-allele-order \ 
  --make-bed \
  --out <out_prefix>
```
expected output files: \
`.bed`: binary file that contains genotype information. \
`.bim`: tab-delimited text file that always accompanies a .bed genotype file. It contains variant information, has no header line, and one line per variant with the following six fields: \
	- Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name \
	- Variant identifier \
	- Position in morgans or centimorgans \
	- Base-pair coordinate (1-based) \
	- Allele 1 (usually minor) \
	- Allele 2 (usually major) \
`.fam`: tab-delimited text file that always accompanies a .bed genotype file. It contains sample information, has no header line, and one line per sample with the following six fields:
	- Family ID ('FID') \
	- Within-family ID ('IID'; cannot be '0') \
	- Within-family ID of father ('0' if father isn't in dataset) \
	- Within-family ID of mother ('0' if mother isn't in dataset) \
	- Sex code ('1' = male, '2' = female, '0' = unknown) \
	- Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control) 

The flag `--make-bed` tells PLINK to output the `*.bed`, `*.bim`, `*.fam` fileset, called the PACKEDPED format. \
There are other PLINK formats but this is the best for working with PLINK and EIGENSOFT downstream. \
Also note that the files correspond to each other so you cannot manually filter one of them without filtering the fileset. \
e.g. You could rename the variant IDs as long as the same number of variants are in the `*.bim` file, but not remove or add variants.

## Subset by individuals in PLINK
Subset bed files to required individuals:
```
module load plink/1.90beta-4.4-21-May

plink --bfile <input_fileset_prefix> \
	--keep-allele-order \
	--allow-no-sex \
	--keep keep_list.txt \
	--make-bed \
	--out <output_fileset_prefix>
```
Where keeplist.txt has one individual per row, the first and second column from the `*.fam` file

## Merge datasets in PLINK

```
plink --bfile <FIRST_input_fileset_prefix> \
	--keep-allele-order \
	--allow-no-sex \
	--merge-list mergelist.txt \
	--make-bed \
	--out <output_fileset_prefix>
```
Where mergelist.txt has the format:
```
SECOND_input.bed SECOND_input.bim SECOND_input.fam
THIRD_input.bed THIRD_input.bim THIRD_input.fam
```

## Convert PLINK to Eigenstrat format
Use EIGENSOFT's CONVERTF for converting formats. \
CONVERTF manual: https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README \
The syntax to use convertf is `convertf -p parfile`
Where the parfile should be named `par.PACKEDPED.EIGENSTRAT.<name>` \
With the following format:
```
genotypename:    <in>.bed
snpname:         <in>.bim
indivname:       <in>.fam
outputformat:    EIGENSTRAT
genotypeoutname: <out>.geno
snpoutname:      <out>.snp
indivoutname:    <out>.ind
```

## Subset by individuals in EIGENSTRAT
Edit the `*.ind` file and change the population to "ignore" for samples you want to remove from the dataset. \
Then use convertf to convert EIGENSTRAT to EIGENSTRAT format and the output will contain your subsetted individuals. \
The parfile will have the name `par.EIGENSTRAT.EIGENSTRAT.<name>`

## Merge datasets in EIGENSTRAT
Use mergeit, syntax is `mergeit -p parfile`. \
mergeit documentation: https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README \
`*.parfile` format:
```
geno1: <input1>.geno
snp1:  <input1>.snp
ind1:  <input1>.ind
geno2: <input2>.geno
snp2:  <input2>.snp
ind2:  <input2>.ind
genooutfilename: <output>.geno
snpoutfilename:	<output>.snp
indoutfilename:	<output>.ind
outputformat:	EIGENSTRAT
docheck:	YES
hashcheck:	YES
```
NB** in the official mergeit documentation, this parfile is incorrect. \
The documentation reads `genotypeoutname` `snpoutname` `indivoutname`, instead of what is in the above example. \

## Running ADMIXTURE
Best to pseudohaploidise data if low-coverage or ancient. \
The ADMIXTURE manual says a minumum of 10,000 markers are necessary for comparing populations betwene continents, \
but at least 100,000 are better for comparing within a continent. \
Requires the `*.bed`, `*.bim`, `*.fam` fileset in the working directory, and then the `*.bed` file is called in the script \
```
ml Admixture/1.3.0

cd <path_to_output_directory>

BED_FILE=<path>/<input_file>.bed

for K in {3..12}; do
		admixture -C 100 -j2 -s time --cv $BED_FILE $K | tee admixture_5_log${K}.out \
done
```
`-C 100` stops the algorithm after 100 iterations \
`-j2` specifies 2 threads \
`-s time` generates a random seed based on the time \
`--cv` Means cross-validation will be calculated \

ADMIXTURE 1.3.0 manual: https://vcru.wisc.edu/simonlab/bioinformatics/programs/admixture/admixture-manual.pdf \
Useful ADMIXTURE tutorial: https://gaworkshop.readthedocs.io/en/latest/contents/07_admixture/admixture.html \

## Plotting ADMIXTURE in R
Download the *.Q file for each K value generated \
Depending on the populations in your `*.ind` file, you may want to download it as well to use as population labelling, \
or write your own file, as long as the order of populations corresponds to the `*.ind` file, cbind will work to label the right sample. \
R code to plot admixture for K=6:
```
#read in data
tbl <- read.table("<path_to_file>/<name>.6.Q")

#read in table for labelling
labTable <- read.table("<path_to_file>/*.labels.txt", col.names=c("Code", "ID", "XLabel"))

#this binds (column-wise) your admixture proportions and your population names taken from the .labels.txt file
mergedAdmixtureTable = cbind(tbl, labTable)

#this orders your table by your populations
ordered <- mergedAdmixtureTable[order(mergedAdmixtureTable$XLabel),]

#plot
barplot(t(as.matrix(subset(ordered, select=V1:V6))), 
        col=rainbow(6), space =0.05, border=NA, ylab="Ancestry", xlab="K=6", names.arg=ordered$XLabel, las = 3)
```
the cv (Cross-Validation) value for each value of K will be buried in the `*log*.out` file \
usually `tail *out` will print the sections you need to stdout

## Running SmartPCA
Useful Tutorial on running PCA: https://gaworkshop.readthedocs.io/en/latest/contents/05_pca/pca.html#running-smartpca \
Syntax to run smartPCA is `smartpca -p parfile`. \
SmartPCA documentation: https://github.com/argriffing/eigensoft/blob/master/POPGEN/README \
Inputs: EIGENSTRAT fileset, and a `*.poplist.txt` file containing one population per line, of the populations in the `*.ind` file. \
Only populations used to calculate eigenvectors should be in this file, then other populations will be automatically projected. \

Format of the parfile, named `*.smartpca.params.txt`. 
```
genotypename:	<path>/<input_name>.geno
snpname:	<path>/<input_name>.snp
indivname:	<path>/<input_name>.ind
evecoutname:	<path>/<ouput_name>.pca.evec.txt
evaloutname:	<path>/<ouput_name>.pca.eval.txt
poplistname:	<path>/<ouput_name>.poplist.txt
lsproject:	YES
outliermode:	2
shrinkmode:	YES
numoutevec:	4
```
`poplistname`: If wishing to infer eigenvectors using only individuals from a 
  subset of populations, and then project individuals from all populations 
  onto those eigenvectors, this input file contains a list of population names,
  one population name per line, which will be used to infer eigenvectors.  It is assumed that the population of each individual is specified in the 
  `*.ind` file.  Default is to use individuals from all populations. \
`outliermode: 2` Do not remove outliers

## Plotting PCA in R
Download the `*evec.txt` and `*ind.` file. \
Depending on your dataset, you will probably want to write a file with the population groups you want to label in the plot, corresponding to the Sample name in the `*ind` and `*evec.txt` files. \
Make sure the order of individuals in all files is the same, so `cbind` works to merge your data. \
R script:
```
#read in data
fn = "<path>/<name>.pca.evec.txt"
evecDat = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", 
			"PC7", "PC8", "PC9", "PC10", "Population"))

#read in population groups to table for labelling purposes
populations = read.table("<path>/<name>.popLabels.txt",col.names=c("Sample", "Label"))
ind = read.table("<path>/<name>.ind", 
                col.names=c("Sample", "Sex", "Population"))

#bind population colums (horizontally) for labelling first
mergedPopDat = cbind(ind, populations)
#merge with evec data
mergedEvecDat3 = merge(mergedPopDat3, evecDat3, by="Sample")

#plot
plot(mergedEvecDat$PC1, mergedEvecDat$PC2, col=mergedEvecDat$Label, 
			pch=as.integer(mergedEvecDat$Label) %% 24, xlab="PC1", ylab="PC2")
legend("topright", xpd=TRUE, legend=levels(mergedEvecDat$Label), 
			col=1:length(levels(mergedEvecDat$Label)), pch=1:length(levels(mergedEvecDat$Label)))
```

## Miscellaneous Useful commands

Renaming SNP ID from the rsID to "CHR_SITE" \
In `*.bim` files:
```
awk '{print $1, "\t", $1"_"$4, "\t", $3, "\t", $4, $5, "\t", $6}' <old>.bim > <new>.2.bim
```
In `*.snp` files:
```
awk '{print $2"_"$4, "\t", $2, "\t", $3, "\t", $4, $5, $6}' <old>.snp > <new>.snp
```

Removing rows in a text file by duplicates in a specified colums.\
e.g. to remove rows with duplicats in column 2:
```
awk '!seen[$2]++' in.txt > out.txt
```

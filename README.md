# Analysis-Scripts
Scripts and notes on how to analyse processed NGS data

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
S
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

## Merge datasets in EIGENSTRAT



## Runing Admixture
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

## Plotting admixture in R
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

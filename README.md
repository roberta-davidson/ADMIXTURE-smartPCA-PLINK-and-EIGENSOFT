# ADMIXTURE-PCA-PLINK-EIGENSTRAT
- Scripts and notes on how to work with PLINK files, EIGENSTRAT files, converting between them. \
- How to run and plot ADMIXTURE and smartPCA. 

## Using PLINK
PLINK is in general very annoying, reccommend to manipulate data in VCF or in EIGENSTRAT formats where possible. \
There are many functions PLINK will do to your data by default, so find the flags necessary to turn off these functions. \
Some useful ones I use: \
- `--keep-allele-order`	Use this EVERY SINGLE TIME you call a plink command, otherwise the order of Allele1 and Allele2 may (or probably will) be flipped in your data. \
- `--allow-no-sex` 	PLINK will default to removing individuals that have unassigned sex, use this to force it to keep them. \
- `--snps-only` 		Removes indels from your variant data and keeps only snps \
- `--biallelic-only`	Removes sites with 2+ alleles \
- `--indiv-sort 0` PLINK default re-orders your data by individual name, this keeps them the same order as the `*.fam` file \
- `--geno xx` removes sites with missingness greater than a given thrshold. PLINK by default filters snps with >0.1 missingness, so use `--geno 1.0` to keep all sites. `--geno 0.999999`	Removes sites no data \
- `--mind` similar to geno, but sets a threshold of missingness per individual. \
- `--extract`/`--exclude` Extracts or exlcludes variants based on a .txt file list of all variant IDs
- `--keep`/`--remove` keep or remove individuals based on a supplied list in a .txt file with corresponding family ID nd within family IDs (or population & individual names). 

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
	- Variant identifier (either rsID or "CHR_POS")\
	- Position in morgans or centimorgans (can leave as 0.0)\
	- Base-pair coordinate (1-based) \
	- Allele 1 (usually minor) \
	- Allele 2 (usually major) \
Example:
```
1	1_752566     	0.0      	752566 G A
1	rs12124819     	0.020242       	776546 A G
8	8_129184555     1.291846    	129184555 C T
```
`.fam`: tab-delimited text file that always accompanies a .bed genotype file. It contains sample information, has no header line, and one line per sample with the following six fields:
	- Family ID ('FID') \
	- Within-family ID ('IID'; cannot be '0') \
	- Within-family ID of father ('0' if father isn't in dataset) \
	- Within-family ID of mother ('0' if mother isn't in dataset) \
	- Sex code ('1' = male, '2' = female, '0' = unknown) \
	- Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control) 
Example:
```
CentralCoast 	I0971 0 0 2 1
SouthCentralHighlands 	PML5 0 0 0 1
NorthChile 	I2540 0 0 1 1
```

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
Use poplistname option in convertf \
Then use convertf to convert EIGENSTRAT to EIGENSTRAT format and the output will contain your subsetted individuals. \
The syntax to use convertf is `convertf -p parfile`
The parfile will have the name `par.EIGENSTRAT.EIGENSTRAT.<name>` 
Example:
```
genotypename:    <input>.geno
snpname:         <input>.snp
indivname:       <input>.ind
outputformat:    EIGENSTRAT
genotypeoutname: <subset>.geno
snpoutname:      <subset>.snp
indivoutname:	 <subset>.ind
poplistname:	 poplist_<subset>.txt
```
Where the file you give to poplistname has been written to include populations (1 per line) from the `.ind` file that you want to extract. \
Depending on your dataset this may be as simple as listing some popualtions, or you may have to set the population column of the individuals you are removing to 'ignore', and edit the poplist accordingly. 

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
Also good practice to prune for linkage disequilibrium \
```
plink --bfile <data> --keep-allele-order --allow-no-sex --indep-pairwise 50 10 0.1 --make-bed --out <data_LDpruned>
```

The ADMIXTURE manual says a minumum of 10,000 markers are necessary for comparing populations betwene continents, \
but at least 100,000 are better for comparing within a continent. \
Requires the `*.bed`, `*.bim`, `*.fam` fileset in the working directory, and then the `*.bed` file is called in the script \
```
ml Admixture/1.3.0

cd <path_to_output_directory>

BED_FILE=<path>/<input_file>.bed

for K in {2..12}; do
		admixture -C 100 -j2 -s time --cv $BED_FILE $K | tee admixture_log${K}.out
done
```
`-C 100` stops the algorithm after 100 iterations \
`-j2` specifies 2 threads \
`-s time` generates a random seed based on the time \
`--cv` Means cross-validation will be calculated \

ADMIXTURE 1.3.0 manual: https://vcru.wisc.edu/simonlab/bioinformatics/programs/admixture/admixture-manual.pdf \
Useful ADMIXTURE tutorial: https://gaworkshop.readthedocs.io/en/latest/contents/07_admixture/admixture.html \

## Supervised ADMIXTURE
If something is known about the relationship between populations, can select populations as fixed source groups of ancestry and infer those ancestries in you test individuals. \
Requires the flag --supervised and an additional `*.pop` file (with matching prefix), this has the same number of lines as the `*.fam`, one line per individual with the population names if they denote a fixed ancestry and those you wish to infer ancestry for have a "-" instead. e.g.
```
Fam1
Fam2
-
-
-
```
Where Fam1 and Fam2 are fixed ancestries and ancestral proportion will be inferred for Ind3-Ind5. \
Note ADMIXTURE will only run for K=(number of fixed ancestries), it will not infer more ancestries at higher K values than the number you specify in the `*.pop` file \
Script:
```
ml Admixture/1.3.0

cd <path_to_output_directory>

BED_FILE=<path>/<input_file>.bed

for K in {3..12}; do
		admixture -C 100 -j2 -s time --supervised --cv $BED_FILE $K | tee admixture_5_log${K}.out \
done
```

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

## Plotting Admixture with PONG
Install: 
```
pip install pong
```
PONG Manual: http://brown.edu/Research/Ramachandran_Lab/files/pong/pong-manual.pdf

Pong command:
```
/Users/robertadavidson/miniconda3/envs/pong/bin/pong -m filemap.txt -n pop_order.txt -i ind2pop.txt
```

Write a filemap text file that gives paths to the matrix files that are the output of Admixture and similar softwares. \

`Column 1` 
The runID, a unique label for the Q matrix (e.g. the string “run5_K7”). Note: A runID must begin with a letter (A-Z/a-z), followed by any number of hyphens (-),
underscores (_), letters, or numbers. Other characters are not allowed in runIDs. Hashmarks (#) can be used in the filemap to indicate the start of a comment.

`Column 2`
The K value for the Q matrix. Each value of K between Kmin and Kmax must be represented by at least one Q matrix in the filemap; if not, pong will abort.

`Column 3`
The path to the Q matrix, relative to the location of the filemap. Thus, if the filemap is in the same directory as the Q matrix file, this is just the name of the Q matrix file. Note that the metadata provided in the filemap allow the user to apply pong to Q matrices in multiple directories in the user’s computer. The path cannot contain a hashmark (#) because it will be interpreted as a comment. \
Example
```
k2r1	2	data/run1/pruned_filtered_1kg_phase3.2.Q
k2r2	2	data/run2/pruned_filtered_1kg_phase3.2.Q
k2r3	2	data/run3/pruned_filtered_1kg_phase3.2.Q
k2r4	2	data/run4/pruned_filtered_1kg_phase3.2.Q
k2r5	2	data/run5/pruned_filtered_1kg_phase3.2.Q
k2r6	2	data/run6/pruned_filtered_1kg_phase3.2.Q
k2r7	2	data/run7/pruned_filtered_1kg_phase3.2.Q
k2r8	2	data/run8/pruned_filtered_1kg_phase3.2.Q
k3r1	3	data/run1/pruned_filtered_1kg_phase3.3.Q
k3r2	3	data/run2/pruned_filtered_1kg_phase3.3.Q
k3r3	3	data/run3/pruned_filtered_1kg_phase3.3.Q
k3r4	3	data/run4/pruned_filtered_1kg_phase3.3.Q
k3r5	3	data/run5/pruned_filtered_1kg_phase3.3.Q
k3r6	3	data/run6/pruned_filtered_1kg_phase3.3.Q
k3r7	3	data/run7/pruned_filtered_1kg_phase3.3.Q
k3r8	3	data/run8/pruned_filtered_1kg_phase3.3.Q
k4r1	4	data/run1/pruned_filtered_1kg_phase3.4.Q
k4r2	4	data/run2/pruned_filtered_1kg_phase3.4.Q
k4r3	4	data/run3/pruned_filtered_1kg_phase3.4.Q
k4r4	4	data/run4/pruned_filtered_1kg_phase3.4.Q
k4r5	4	data/run5/pruned_filtered_1kg_phase3.4.Q
```
Write an ind2pop.txt file that has one population label per line, corresponding to the order of samples in the Admixture output matrices.
Example:
```
Mbuti
Japanese
French
French
Mayan
Japanese
Mayan
```

optional: Write a file to determine the order populations are plotted in the PONG output bar graph. \
`Column 1` labels correspond to the ind2pop file
`Column 2` (Optional) Longer labels for the population
Example
```
Mbuti		Label 1
French		Label 2
Japanese	Label 3
Mayan		Label 4
```
Typical PONG output:
<img width="1019" alt="image" src="https://user-images.githubusercontent.com/78726635/156501053-c0208021-5b12-418c-91aa-b5c597adc7c9.png">

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
numoutevec:	10
```
`poplistname`: If wishing to infer eigenvectors using only individuals from a 
  subset of populations, and then project individuals from all populations 
  onto those eigenvectors, this input file contains a list of population names,
  one population name per line, which will be used to infer eigenvectors.  It is assumed that the population of each individual is specified in the 
  `*.ind` file.  Default is to use individuals from all populations. \
`outliermode: 2` Do not remove outliers

## Running emu for PCA
paper: https://academic.oup.com/bioinformatics/article/37/13/1868/6103565?login=true \
EMU Github & installation instructions: https://github.com/Rosemeis/emu \
EMU works with PLINK binary filesets and has inbuilt maf filtering. \
It better compensates for missingness typical of aDNA \
Requires python 2.7 so create conda environment for that: `conda create -n "emu" python=2.7` \
To install & build emu:
```
git clone https://github.com/Rosemeis/emu.git
cd emu
python setup.py build_ext --inplace
pip3 install -e .
```

Script to run emu on dataset:
```
conda activate emu

DATA=$1

emu --plink ${DATA} \
	--n_eig 10 \
	--threads 8 \
	--maf 0.01 \
	--out ${DATA} \
	--maf_save \
	--sites_save
```


## Plotting PCA in R
Download the `*evec.txt` and `*ind.` file. \
Depending on your dataset, you will probably want to write a file with the population groups you want to label in the plot, corresponding to the Sample name in the `*ind` and `*evec.txt` files. \
Make sure the order of individuals in all files is the same, so `cbind` works to merge your data. \
R scripts in base R and ggplot (You choose):
base R plot:
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
ggplot with faceting:
```
library(ggplot2)
setwd("/Users/robbi/Box/Robbi_PhD/02_Inka_Royal_Ancestry/IncaModern/EAG/") #set working directory

fn = "/Users/robbi/Box/Robbi_PhD/02_Inka_Royal_Ancestry/IncaModern/EAG/pca_1.pca.evec.txt"
evecDat1 = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "Pop"))

#make colour palette
colours <- c("blue", "#E2AC36", "#B3493A", "#47271C", "#BCAB7E","#BCAB7E")
colours <- colours[as.factor(evecDat1$Pop)]  

#define symbols
shapes <- c(0,0,3,0,0)
shapes <- shapes[as.factor(evecDat1$Pop)]  

#base plot 1v2
plot_1v2 <- (
  plot(evecDat1$PC1, evecDat1$PC2, col = colours, xlab="PC1", ylab="PC2", pch = shapes, cex = 2, lwd=1.5,
       cex.axis=1.3, cex.lab=1.2, las=1))
legend("topleft", legend = levels(as.factor(evecDat1$Pop)), col = c("blue", "#E2AC36", "#B3493A", "#47271C", "#BCAB7E","#BCAB7E"),
       pch = c(0,0,3,0,0), cex = 1.2, bty="o", pt.cex=2, pt.lwd=1.5, scale_fill_discrete(labels = labels))
#ggplot 1v2
PC1_2 <- ggplot(evecDat1, #dataset to plot
              aes(x = PC1, #x-axis is petal width
                  y = PC2, #y-axis is petal length
                  color = Pop)) + #each species is represented by a different shape
  geom_point() + #default scatter plot
  theme_light() + #light theme with no grey background and with grey lines and axes
  scale_x_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on top 
  scale_y_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on right
  theme(panel.grid.major = element_blank(), #remove major gridline
        panel.grid.minor = element_blank(), #remove minor gridline
#        legend.justification = c(0, 0), #justification is centered by default, c(1, 0) means bottom right
#        legend.position = c(0.97, 0.01), #position relative to justification
#        legend.background = element_rect(color = "grey"), #legend box with grey lines
#        legend.text = element_text(face = "italic"), #since the legend is species names, display in italics
        axis.title.x.top = element_blank(), #remove top x-axis title
        axis.text.x.top = element_blank(), #remove labels on top x-axis
        axis.title.y.right = element_blank(), #remove right y-axis title
        axis.text.y.right = element_blank()) + #remove labels on right y-axis
  scale_shape(labels = c("Mbuti", "Karitiana", "Han", "Basque", "IncaDescendant"), #edit legend labels
              solid = FALSE) + #force hollow points to increase clarity in case of overlaps
  labs(x = "PC1 (7.633 %)", #edit x-axis label
       y = "PC2 (3.975 %)", #edit y-axis label
       shape = "Population") #edit legend title
PC1_2 #show plot

#ggplot 1v2
PC1_3 <- ggplot(evecDat1, #dataset to plot
                aes(x = PC1, #x-axis is petal width
                    y = PC3, #y-axis is petal length
                    color = Pop)) + #each species is represented by a different shape
  geom_point() + #default scatter plot
  theme_light() + #light theme with no grey background and with grey lines and axes
  scale_x_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on top 
  scale_y_continuous(sec.axis = dup_axis()) + #add duplicated secondary axis on right
  theme(panel.grid.major = element_blank(), #remove major gridline
        panel.grid.minor = element_blank(), #remove minor gridline
        legend.justification = c(0, 0), #justification is centered by default, c(1, 0) means bottom right
        #        legend.position = c(0.97, 0.01), #position relative to justification
        legend.background = element_rect(color = "grey"), #legend box with grey lines
        #        legend.text = element_text(face = "italic"), #since the legend is species names, display in italics
        axis.title.x.top = element_blank(), #remove top x-axis title
        axis.text.x.top = element_blank(), #remove labels on top x-axis
        axis.title.y.right = element_blank(), #remove right y-axis title
        axis.text.y.right = element_blank()) + #remove labels on right y-axis
  scale_shape(labels = c("Mbuti", "Karitiana", "Han", "Basque", "IncaDescendant"), #edit legend labels
              solid = FALSE) + #force hollow points to increase clarity in case of overlaps
  labs(x = "PC1 (7.633 %)", #edit x-axis label
       y = "PC3 (2.650 %)", #edit y-axis label
       shape = "Population") #edit legend title
PC1_3 #show plot

#composite 
library(patchwork)
PC1_2_patch <- PC1_2 + #prepare a version of gg3 for patchwork design
      theme(legend.position="none") #remove legend
PC1_2_patch #show plot

patch <- PC1_2_patch + PC1_3 + #assemble patchwork
  plot_layout(widths = c(3, 3)) #set width ratio of the 3 panels
#  plot_annotation(tag_levels = "A") #add plot labels (uppercase Latin letters)
patch #show plot

ggsave("patch.pdf", width = 12, height = 6) #save in pdf format with size 12 x 6 in
```
## Miscellaneous Useful commands

Renaming SNP ID from the rsID to "CHR_SITE" \
In `*.bim` files:
```
awk '{print $1,$1"_"$4,$3,$4,$5,$6}' <old>.bim > <new>.2.bim
```
In `*.snp` files:
```
awk '{print $2"_"$4,$2,$3,$4,$5,$6}' <old>.snp > <new>.snp
```
Then to replace spaces with tabs:
```
sed -i 's/ /\t/g' <file>
```

Removing rows in a text file by duplicates in a specified colums.\
e.g. to remove rows with duplicats in column 2:
```
awk '!seen[$2]++' in.txt > out.txt
```
Editing `.ind` file to set population name to 'ignore' for individuals other than ones you want to keep. \
(Only really practical if subsetting for a small number of individuals)
```
awk '{if ($1=="Sample1"||$1=="Sample2"||$1=="Sample3") print $0; else print $1, $2, "ignore"}' v44.3_1240K_public.ind > v44.3_1240K_public.subset.ind

```
Find and replace strings in text file. `\b` denotes word boundary
```
gsed -i 's/\b<OLD_STRING>\b/<NEW_STRING>/g' <file>.txt
```

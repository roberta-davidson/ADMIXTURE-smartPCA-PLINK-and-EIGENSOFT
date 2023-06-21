ml plink/1.90beta-4.4-21-May

in1=$1
in2=$2
out=$3

plink --bfile ${in1} #input fileset
  --bmerge ${in2} #filest to merge with
  --merge-mode 5 # "always overwrite" mode to resolve conflicts during merging - check documentation for others
  --keep-allele-order # prevent plink from re-ordering alleles for it's own evil fun
  --allow-no-sex #block error about individuals with no sex
  --mind 1.0 #prevent automatc removal of inds with low SNPs
  --geno 1.0 #prevent automatc removal of SNPs with high missingness
  --indiv-sort 0 # prevent automatic .ind re-ordering
  --make-bed --out ${out} #output fileset

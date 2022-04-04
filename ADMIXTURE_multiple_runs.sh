#load modules
ml Admixture/1.3.0
ml plink

# Define variables
DATA=$1

#make major analysis directory
mkdir ./Admixture_${DATA}
cd ./Admixture_${DATA}

#first filter for missing snps sites
plink --bfile ../${DATA} --allow-no-sex --keep-allele-order --geno 0.99999 --mind 1.0 --make-bed --out ./${DATA}_2
rm *nosex
echo "#### plink filtering done ####"

#make directories for each run and copy plink dataset and admixture script into each
for R in {1..10}; do
mkdir run_${R}
cp ${DATA}_2 run_${R}
cp Admixture.sh run_${R}
done
echo '#### directories setup ####'

#Run Admixture 
for R in {1..10}; do
cd run_${R}
sbatch Admixture.sh ${DATA}_2
done

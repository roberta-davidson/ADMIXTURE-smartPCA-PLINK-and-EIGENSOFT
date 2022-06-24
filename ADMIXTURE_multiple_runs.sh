#load modules
ml Admixture/1.3.0
ml plink

# Define variables
DATA=$1
RUNS=10

#make major analysis directory
mkdir ./Admixture_${DATA}
cd ./Admixture_${DATA}

#first filter for missing snps sites
plink --bfile ../${DATA} --allow-no-sex --keep-allele-order --geno 0.99999 --mind 1.0 --make-bed --out ./${DATA}_2
rm *nosex
echo "#### plink filtering done ####"

#### ADD LD pruning step####

#make directories for each run and copy plink dataset and admixture script into each
for R in {1..10}; do
	mkdir run_${R}
	cp ${DATA}_2.* run_${R}
	cp ../Admixture.sh run_${R}
done
echo '#### directories setup ####'

#Run Admixture in parallel for each independent run
parallel=$RUNS
for R in {1..10}; do
	((i=i%parallel)); ((i++==0))
	cd run_${R}
	sbatch Admixture.sh ${DATA}_2 
	cd .. &
done

# Write cross validation text file for all runs
echo "writing CV file for each run"

for R in {1..10}; do
	cd run_${R}
	grep "CV" admixture_data19_2_K_* > Run${R}_CV.txt
	sed -i 's/:/\t/g' Run${R}_CV.txt
	sed -i 's/=/\t/g' Run${R}_CV.txt 
	sed -i 's/)//g' Run${R}_CV.txt 
	#print columns
	awk -v OFS='\t' '{print $5,$6}' Run${R}_CV.txt > Run${R}_CV_2.txt
	#sort
	sort -n Run${R}_CV_2.txt > Run${R}_CV_sort.txt
	#make header
	echo -e "K\trun_${R}" > Run${R}_CV_header.txt 
	#CV plus header
	cat Run${R}_CV_header.txt Run${R}_CV_sort.txt > Run${R}_CV.txt
	#cleanup
	rm Run${R}_CV_sort.txt Run${R}_CV_2.txt Run${R}_CV_header.txt
	cd ..
done
echo "done"

#write master CV file
echo "writing master CV file"
paste run_*/Run*_CV.txt > Master_CV.txt
sed -i 's/ /\t/g' Master_CV.txt
echo "done"

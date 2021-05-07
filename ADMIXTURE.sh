ml Admixture/1.3.0

cd <path_to_output_directory>

BED_FILE=<path>/<name>.bed

for K in {3..12}; do
	admixture -C 100 -j2 -s time --cv $BED_FILE $K | tee <name>_log${K}.out
done

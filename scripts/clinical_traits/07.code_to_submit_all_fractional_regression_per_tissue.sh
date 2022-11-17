#!/bin/bash
#This code uses as input the tissue name
#It executes voom limma without interaction and with interactions per each of the diseases


folder=/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/03_Models/Tissues/$1
for file in "$folder"/*; do
if [[ $file == *"SampleMetadata"* ]]; then
	acronyms="${file#*.}"
	acronyms=${acronyms#*.}
	acronyms=${acronyms%.*}
	echo $acronyms
	echo $1
	sbatch --constraint=highmem --cpus-per-task=16 07.to_submit_fractional_regression_per_tissue.sh $acronyms $1 
fi
done


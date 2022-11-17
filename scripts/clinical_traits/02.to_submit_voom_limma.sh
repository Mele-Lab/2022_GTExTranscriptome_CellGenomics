#!/bin/sh

echo $1
cat /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissues.txt | while read line; do Rscript /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/03_Models/02.voom_limma.Disease_publication.R -t $line -d $1 -o /gpfs/projects/bsc83/Projects/GTEx_v8/Jose/01_Overview/Final_Diseases/; done; 

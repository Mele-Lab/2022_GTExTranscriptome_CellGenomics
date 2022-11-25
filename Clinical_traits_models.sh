
#Preprocessing histological data, where ${sample_metadata.csv} is the file publicly available downloaded from gtex.

Rscript scripts/00.parsing_histological_data.R -f ${sample_metadata.csv} #histological_data.csv file is created

#Metadata per donor disease
Rscript scripts/01.Donor_disease_phenotypes.R -a /gpfs/projects/bsc83/Data/GTEx_v8/Annotations/ -t ../00_Data/Tissue_info.rds -d /gpfs/projects/bsc83/Projects/GTEx_v8/ #It uses the protected file obtained from dbGAP

#Metadata per sample disease:
Rscript scripts/01.Sample_disease_phenotypes.R -f ../00_Data/histological_data.csv -t ../00_Data/Tissue_info.rds -d /gpfs/projects/bsc83/Projects/GTEx_v8/ 

#Merging the names of the sample and donor diseases (using the byproducts of the two previous scripts)

cat sample_diseases.txt > diseases.txt
cat donor_diseases.txt >> diseases.txt


#DEA submitted to a cluster with SLURM. We could run an array instead
cat diseases.txt | while read line;
do sbatch --time=00:45:00 scripts/02.to_submit_voom_limma.sh $line; done #It calls scripts//02.voom_limma.Disease_publication.R

Rscript scripts/02.Number_of_DEG_per_tissue.R -d /gpfs/projects/bsc83/Projects/GTEx_v8/ #It outputs heatmap_data_PEER.RData

#Final models
Rscript scripts/03.Tissue_names.R #It creates Tissues.txt

tissue_file="/Documents/mn4/Jose/00_Data/Tissues.txt"
while read tissue; do # reading each line
Rscript scripts/04.PrepareMetadata_per_tissue.R $tissue
folder=~/Documents/mn4/Jose/03_Models/Tissues/$tissue
for file in "$folder"/*; do
if [[ $file == *"SampleMetadata"* ]]; then
	acronyms="${file#*.}"
	acronyms=${acronyms#*.}
	acronyms=${acronyms%.*}
	echo ${acronyms}
	Rscript scripts/02.voom_limma.Disease_publication.R -t $tissue -d "$acronyms"
fi
done
done < $tissue_file


#Hier part
Rscript scripts/05.Residuals_per_tissue.R #It may take long, we recommend a cluster environment with --constraint=highmem --cpus-per-task=48

#Saving DEA results into the table available in zenodo
Rscript scripts/06.differential_expression_analysis_tables.R

#Splicing. The code runs DSA and hier part
ls Tissues/ | while read line;
do scripts/07.code_to_submit_all_fractional_regression_per_tissue.sh $line; done

#Saving splicing results into the table available at zenodo
Rscript scripts/08.differential_splicing_analysis.R


#xCell
Rscript scripts/09.xCell.Sex_paper_model.Only_2_covariates_per_tissue.R

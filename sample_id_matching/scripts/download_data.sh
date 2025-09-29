#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --mem=4G
#SBATCH --time=1-12:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=sra-toolkit
#SBATCH --output=download_data_%J.log
#SBATCH --error=download_data_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

# options: central_nervous_system_tumor, hgsoc_ariel, high_grade_glioma, low_grade_glioma, wilms_tumor
dataset=$1
# options: bulk, single-cell, single-nucleus, all
data_type=$2

module load sra-toolkit

mkdir -p "../data/${dataset}/${data_type}"
download_location="/scratch/alpine/$USER/hgsoc/simulations/data/${dataset}/${data_type}"
final_fastq_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/${dataset}/${data_type}"

# Use .txt file of accession numbers to prefetch .sra files
# prefetch \
#    --option-file "../data/accession_lists/SRR_Acc_List_${data_type}_${dataset}.txt" \
#    --output-directory "${download_location}"

# Download .fastq files from .sra files using fasterq-dump
# for accession in $(cat "../data/accession_lists/SRR_Acc_List_${data_type}_${dataset}.txt")
# do
#     if [ ! -f "${final_fastq_location}/${accession}_*.fastq.gz" ]
#     then 
#         fasterq-dump \
#             --split-files \
#             --outdir "${download_location}" \
#             "${accession}"
#     fi
# done

for accession in $(cat "../data/accession_lists/SRR_Acc_List_${data_type}_${dataset}.txt")
do
    if ls "${final_fastq_location}/${accession}"*".fastq.gz" >/dev/null 2>&1
    then
        echo "${accession} already downloaded"
    else
        # Force re-download if files are missing
        rm -r "${download_location}/${accession}"

        prefetch \
            "${accession}" \
            --output-directory "${download_location}"

        fasterq-dump \
            --split-files \
            --outdir "${download_location}" \
            "${accession}"
    fi
done

# Zip .fastq files
for file in "${download_location}/"*".fastq"
do 
    gzip "$file"
done

# Move .fastq files to final location
mv "${download_location}/"*".fastq.gz" "${final_fastq_location}/"
# Clean up .sra files to save space
# rm -r "${download_location}/"*".sra"

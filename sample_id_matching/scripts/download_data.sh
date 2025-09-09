#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=sra-toolkit
#SBATCH --output=download_data_%J.log
#SBATCH --error=download_data_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

# options: sn_wilms_tumor, sn_cnst, sc_high_grade_glioma, low_grade_glioma, hgsoc, 
# bulk_wilms_tumor, bulk_high_grade_glioma, bulk_cnst
dataset=$1

module load sra-toolkit

mkdir -p "../data/${dataset}"
download_location="/scratch/alpine/$USER/hgsoc/simulations/data/${dataset}"
final_fastq_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/${dataset}"

# Use .txt file of accession numbers to prefetch .sra files
# prefetch \
#     --option-file "../data/SRR_Acc_List_${dataset}.txt" \
#     --output-directory "${download_location}"

# Download .fastq files from .sra files using fasterq-dump
fasterq-dump \
    --split-files \
    --threads 3 \
    --outdir "${download_location}" \
    "${download_location}/SRR"*"/"*".sra"

# Compress .fastq files to save space
for file in "${download_location}/"*".fastq"
do 
    gzip "$file"
done

# Move .fastq files to final location
# mv "${download_location}/"*".fastq.gz" "${final_fastq_location}/"
# Clean up SRR folders with .sra files to save space
# rm -r "${download_location}/SRR"*""

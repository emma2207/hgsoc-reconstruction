#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools
#SBATCH --output=02_bcftools_%J.log
#SBATCH --error=02_bcftools_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load bcftools/1.16

# Script location is /scratch/alpine/$USER/hgsoc/pooled
original_location=`pwd`
data_location="$original_location/bcftools/all"

input_file_bulk="$data_location/bcftools_remove_samples1_bulk_all.vcf"
input_file_diss_bulk="$data_location/bcftools_remove_samples1_diss_bulk_all.vcf"
output_file="$data_location/bcftools_all_samples.vcf.gz"
output_file_bulk="$data_location/bcftools_bulk_filter_mean.vcf"
output_file_diss_bulk="$data_location/bcftools_diss_bulk_filter_mean.vcf" 


# bcftools view -I $input_file_bulk -Oz -o $input_file_bulk.gz
# bcftools view -I $input_file_diss_bulk -Oz -o $input_file_diss_bulk.gz

# bcftools index $input_file_bulk.gz
# bcftools index $input_file_diss_bulk.gz

# bcftools merge -Oz \
#     -o $output_file \
#     $input_file_bulk.gz $input_file_diss_bulk.gz

bcftools view -i 'MEAN(FORMAT/DP)>30' $input_file_bulk.gz > $output_file_bulk.gz
bcftools view -i 'MEAN(FORMAT/DP)>30' $input_file_diss_bulk.gz > $output_file_diss_bulk.gz

bcftools view $output_file_bulk.gz -Ov -o $output_file_bulk
bcftools view $output_file_diss_bulk.gz -Ov -o $output_file_diss_bulk

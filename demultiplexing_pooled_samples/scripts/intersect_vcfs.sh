#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=10
#SBATCH --account=amc-general
#SBATCH --job-name=intersect
#SBATCH --output=intersect_%J.log
#SBATCH --error=intersect_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load bcftools/1.16

# Script location is /scratch/alpine/$USER/hgsoc/pooled
original_location=`pwd`
genomes_location="/projects/$USER/hgsoc/1000_genomes_vcf"
input_file_bulk="$original_location/bcftools/all/bcftools_remove_samples1_bulk_all.vcf.gz"
input_file_diss_bulk="$original_location/bcftools/all/bcftools_remove_samples1_diss_bulk_all.vcf.gz"

# Go through all 1000 genomes files
for file in $genomes_location/*
do
    if [[ "$file" == *.vcf.gz ]]
    then
        filename=$(basename $file)
        echo $filename

        # Get chromosome number
        tmp=${filename#*.chr}
        tmp=${tmp%.*}
        num=${tmp%.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf}
        echo $num

        # Bulk
        mkdir -p bcftools/intersect/bulk/chr$num
        output_dir_bulk="$original_location/bcftools/intersect/bulk/chr$num"
        bcftools isec -p $output_dir_bulk -Oz $input_file_bulk $file

        # Dissociated bulk
        mkdir -p bcftools/intersect/diss_bulk/chr$num
        output_dir_diss_bulk="$original_location/bcftools/intersect/diss_bulk/chr$num"
        bcftools isec -p $output_dir_diss_bulk -Oz $input_file_diss_bulk $file
    fi
done

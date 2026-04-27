#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=split-bam
#SBATCH --output=split-bam_%A_%a.log
#SBATCH --error=split-bam_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load samtools

export pool=$SLURM_ARRAY_TASK_ID

# Split pooled BAM file into sample-specific BAM files based on the read groups (@RG)
current_location=$(pwd)
pooled_bam_location="/scratch/alpine/$USER/hgsoc/pooled/bam/pool${pool}/pooled.bam"
output_location="/scratch/alpine/$USER/hgsoc/pooled/split_bam/pool${pool}"

samtools split -@ 4 -r "RG:Z:sample1" ${pooled_bam_location} ${output_location}/sample1/sample.bam
samtools split -@ 4 -r "RG:Z:sample2" ${pooled_bam_location} ${output_location}/sample2/sample.bam
samtools split -@ 4 -r "RG:Z:sample3" ${pooled_bam_location} ${output_location}/sample3/sample.bam
samtools split -@ 4 -r "RG:Z:sample4" ${pooled_bam_location} ${output_location}/sample4/sample.bam

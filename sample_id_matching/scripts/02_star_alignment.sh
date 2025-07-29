#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --ntasks=6
#SBATCH --account=amc-general
#SBATCH --job-name=star
#SBATCH --output=01_make_bam_file_%A_%a.log
#SBATCH --error=01_make_bam_file_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load star/2.7.10b
module load samtools/1.16.1

dataset=$1
sample=$2
fastq_R1=$3
fastq_R2=$4

index_location="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A/star"
mkdir -p /scratch/alpine/$USER/hgsoc/simulations/aligned_reads/$dataset/$sample
output_location="/scratch/alpine/$USER/hgsoc/simulations/aligned_reads/$dataset/$sample"

STAR \
	--genomeDir ${index_location} \
	--runThreadN 6 \
	--readFilesIn $fastq_R1 $fastq_R2 \
	--outFileNamePrefix ${output_location}/ \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts

samtools index ${output_location}/Aligned.sortedByCoord.out.bam

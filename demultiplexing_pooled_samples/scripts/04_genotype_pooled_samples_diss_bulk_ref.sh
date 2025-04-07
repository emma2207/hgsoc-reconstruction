#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=5G
#SBATCH --time=03:00:00
#SBATCH --ntasks=10
#SBATCH --account=amc-general
#SBATCH --job-name=cellsnp-lite
#SBATCH --output=03_genotype_pooled_%A_%a.log
#SBATCH --error=03_genotype_pooled_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate cellsnp-lite_install

# Based on https://github.com/greenelab/deconvolution_pilot,
# file scripts/genetic_dmx/05_run_cellsnp.sh

export pool=$SLURM_ARRAY_TASK_ID

# Original location is /scratch/alpine/$USER/hgsoc/pooled
mkdir -p cellSNP_diss/pool${pool}

original_location=$(pwd)
pooled_bam_location="$original_location/bam/pool${pool}"
bulk_vcf_location="$original_location/bcftools/pool${pool}"
output_location="$original_location/cellSNP_diss/pool${pool}"

cellsnp-lite \
	-s $pooled_bam_location/pooled.bam \
	-b barcodes/barcodes_pool${pool}.tsv \
	-O $output_location \
	-R $bulk_vcf_location/bcftools_diss_bulk_pool${pool}.vcf \
	-p 10 \
	--minMAF=0.1 \
	--minCOUNT=20 \
	--gzip

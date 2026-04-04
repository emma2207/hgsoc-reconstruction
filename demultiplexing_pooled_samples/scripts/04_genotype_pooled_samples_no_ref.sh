#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=240G
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=4
#SBATCH --account=amc-general
#SBATCH --job-name=cellsnp-lite
#SBATCH --output=03_genotype_pooled_%A_%a.log
#SBATCH --error=03_genotype_pooled_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate cellsnp-lite_install

dataset="hgsoc-new"
pool=10

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
original_location=$(pwd)
pooled_bam_location="$original_location/results/star/${dataset}/pool${pool}"
output_location="$original_location/results/cellSNP/${dataset}/pool${pool}"
barcodes="${original_location}/barcodes_pool${pool}.tsv"

if [ ! -f "$barcodes" ]; then
    echo "Barcode file not found: $barcodes"
    exit 1
fi

mkdir -p ${output_location}

# Remove -1 suffix from barcodes if present (e.g. from 10x v3 chemistry)
if cat $barcodes | head -n 1 | grep -qE -- '-1$'; then
    echo "Removing -1 suffix from barcodes in $barcodes"
    cat $barcodes | sed 's/-1$//' > $output_location/barcodes.tsv
    barcodes="$output_location/barcodes.tsv"
else
    echo "No -1 suffix found in barcodes, using original file: $barcodes"
fi

cellsnp-lite \
	-s "$pooled_bam_location/Aligned.sortedByCoord.out.bam" \
	-b "$barcodes" \
	-O "$output_location" \
	-p 4 \
	--minMAF=0.1 \
	--minCOUNT=20 \
	--gzip

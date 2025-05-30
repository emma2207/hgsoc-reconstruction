#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=128G
#SBATCH --time=00:45:00
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=star
#SBATCH --output=00_create_STAR_index_%J.log
#SBATCH --error=00_create_STAR_index_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load star/2.7.10b

# Based on https://github.com/greenelab/deconvolution_pilot, file scripts/genetic_dmx/00_create_STAR_index.sh
# We will align to the same reference genome that cellranger uses, as obtained from the 10x website:
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest.
# This should only need to be run once, so it will check if the file exists first.

# Current location is /scratch/alpine/elathouwers@xsede.org/hgsoc/scripts
current_location=`pwd`
index_location="/projects/elathouwers@xsede.org/hgsoc/refdata-gex-GRCh38-2024-A"

# Check if STAR index exists
if [ -f "$index_location/star/SAindex" ]; then
	echo "STAR index exists, skipping."
else
	STAR \
		--runMode genomeGenerate \
		--runThreadN 20 \
		--genomeDir $index_location/star \
		--genomeFastaFiles $index_location/fasta/genome.fa \
		--sjdbGTFfile $index_location/genes/genes.gtf
fi

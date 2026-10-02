#!/bin/sh

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=10
#SBATCH --account=amc-general
#SBATCH --job-name=cellsnp-lite
#SBATCH --output=03_cellsnp_%A_%a.log
#SBATCH --error=03_cellsnp_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate cellsnp-lite_install

export pool="$SLURM_ARRAY_TASK_ID"

mkdir -p "cellSNP_bulk_chunk_ribo_specific_refs/pool${pool}"

pooled_bam="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/aligned_reads/hgsoc-new/pooled_single_cell/pool${pool}/pooled.bam"
barcodes_zip=$(ls /pl/active/cgreene-sc-hgsoc/ariel_sc_HGSOC/pooled_sc/Pool${pool}-GEX-CaseyGreene-*/outs/raw_feature_bc_matrix/barcodes.tsv.gz)
bulk_ref="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/variant_calls/hgsoc-new/bulk_chunk_ribo_pool${pool}_ref.vcf.gz"
output_location="cellSNP_bulk_chunk_ribo_specific_refs/pool${pool}"

cellsnp-lite \
	-s "${pooled_bam}" \
	-b "${barcodes_zip}" \
	-O "${output_location}" \
	-R "${bulk_ref}" \
	-p 10 \
	--minMAF=0.1 \
	--minCOUNT=20 \
	--gzip

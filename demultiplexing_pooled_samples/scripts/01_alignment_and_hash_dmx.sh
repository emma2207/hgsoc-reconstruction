#!/bin/sh

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=16
#SBATCH --account=amc-general
#SBATCH --job-name=cellranger_count
#SBATCH --output=01_cellranger_count_%A_%a.log
#SBATCH --error=01_cellranger_count_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

# module use --append /projects/$USER/lmod-files
# module load cellranger/9.0.1
module load cellranger/7.1.0

export pool=$SLURM_ARRAY_TASK_ID
dataset="hgsoc"
datatype="pooled_sc"

# Location is /scratch/alpine/$USER/hgsoc/demultiplexing

cellranger count \
	--id=pool$pool \
	--transcriptome="/projects/$USER/hgsoc/refdata-gex-GRCh38-2024-A" \
	--fastqs="/pl/active/cgreene-sc-hgsoc/ariel_sc_HGSOC/$datatype/fastq" \
	--sample=Pool$pool-GEX \
	--nosecondary
#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=5G
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=vireo
#SBATCH --output=04_vireo_%A_%a.log
#SBATCH --error=04_vireo_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate vireo_install

export pool=$SLURM_ARRAY_TASK_ID
data_type=$1  # "bulk" or "diss_bulk"

# Original location is /scratch/alpine/$USER/hgsoc/pooled
mkdir -p vireo_${data_type}/pool${pool}

original_location=$(pwd)
bulk_vcf_location="${original_location}/bcftools/pool${pool}"
cellsnp_location="${original_location}/cellSNP_${data_type}/pool${pool}"
output_location="${original_location}/vireo_${data_type}/pool${pool}"

# Adjust the number of samples for pools 9 & 10
if [ "$pool" -lt 9 ]; then
	n_samples=4
elif [ "$pool" -lt 10 ]; then
	n_samples=3
else
	n_samples=2
fi

vireo \
	-c ${cellsnp_location} \
	-N ${n_samples} \
	-o ${output_location} \
	-d ${bulk_vcf_location}/bcftools_${data_type}_pool${pool}_rehead.vcf \
	--randSeed=12

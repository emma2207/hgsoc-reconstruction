#!/bin/sh

#SBATCH --array=5,7-9
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=vireo
#SBATCH --output=04_run_vireo_%A_%a.log
#SBATCH --error=04_run_vireo_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate vireo_install

export pool=$SLURM_ARRAY_TASK_ID
data_type=$1  # "bulk" or "diss_bulk"

# Original location is /scratch/alpine/$USER/hgsoc/pooled
mkdir -p vireo_all_ref/${data_type}/pool${pool}

original_location=$(pwd)
bulk_vcf_location="$original_location/bcftools/all"
cellsnp_location="$original_location/cellSNP_no_ref/pool${pool}"
output_location="$original_location/vireo_all_ref/${data_type}/pool${pool}"

# Adjust the number of samples for pools 9 & 10
if [ "$pool" -lt 9 ]; then
	n_samples=4
elif [ "$pool" -lt 10 ]; then
	n_samples=3
else
	n_samples=2
fi

# Select the appropriate bulk VCF file based on data type
if [ "${data_type}" == "diss_bulk" ]; then 
	bulk_vcf_file=$bulk_vcf_location/bcftools_qual_dp_rm_2507_filter_${data_type}_all.vcf
elif [ "${data_type}" == "bulk" ]; then 
	bulk_vcf_file=$bulk_vcf_location/bcftools_qual_dp_filter_${data_type}_all.vcf
else
	echo "Invalid data type specified. Use 'bulk' or 'diss_bulk'."
	exit 1
fi

# Run vireo
vireo \
	-c $cellsnp_location \
	-N $n_samples \
	-o $output_location \
	-d $bulk_vcf_file \
	--randSeed=12

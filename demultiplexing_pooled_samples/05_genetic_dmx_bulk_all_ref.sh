#!/bin/sh

#SBATCH --array=5,7-9
#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
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
sample_to_exclude=$2  # Sample to exclude from the bulk VCF file

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
mkdir -p results/vireo/${data_type}/pool${pool}

original_location=$(pwd)
bulk_diss_polyA_vcf_location="$original_location/../simulations/scripts/experiment/results/1a_vcf/hgsoc-new/real_data/ncells_null/read_depth_bulk_diss_polyA_0_pooled_single_cell_0"
cellsnp_location="$original_location/results/cellSNP/pool${pool}"
output_location="$original_location/results/vireo/${data_type}/pool${pool}"

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
	bulk_vcf_file=$bulk_diss_polyA_vcf_location/bulk_diss_polyA_modality_variants.vcf.gz
elif [ "${data_type}" == "bulk" ]; then 
	echo "Not yet available for bulk data type. Please use 'diss_bulk' for now."
	exit 1
else
	echo "Invalid data type specified. Use 'bulk' or 'diss_bulk'."
	exit 1
fi

# Exclude low quality samples from the bulk VCF file
# Match sample_to_exclude (which is only a partial match) with the sample names in the VCF file and exclude it
full_sample_to_exclude=$(bcftools query -l "$bulk_vcf_file" | grep -E "^HGSOC-${sample_to_exclude}_" || true)

if [ -z "$full_sample_to_exclude" ]; then
	echo "No VCF sample matched sample_to_exclude=${sample_to_exclude}."
	echo "Expected a sample starting with HGSOC-${sample_to_exclude}_."
	exit 1
fi

match_count=$(printf "%s\n" "$full_sample_to_exclude" | wc -l | tr -d ' ')
if [ "$match_count" -ne 1 ]; then
	echo "Expected exactly 1 VCF sample match for sample_to_exclude=${sample_to_exclude}, found ${match_count}."
	echo "Matches:"
	printf "%s\n" "$full_sample_to_exclude"
	exit 1
fi

filtered_bulk_vcf_file=${bulk_vcf_file%.vcf.gz}_rehead.vcf.gz
bcftools view -s ^${full_sample_to_exclude} -Oz -o "$filtered_bulk_vcf_file" "$bulk_vcf_file"
bcftools index -t -f "$filtered_bulk_vcf_file"

# Run vireo
vireo \
	-c $cellsnp_location \
	-N $n_samples \
	-o $output_location \
	-d $filtered_bulk_vcf_file \
	--randSeed=12

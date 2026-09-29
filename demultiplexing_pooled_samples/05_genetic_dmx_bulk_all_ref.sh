#!/bin/sh

#SBATCH --array=2
#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --mem=32G
#SBATCH --time=01:00:00
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
samples_to_exclude=$2  # Sample to exclude from the bulk VCF file (comma separated)

# Original location is /scratch/alpine/$USER/hgsoc/demultiplexing
mkdir -p results/vireo/${data_type}/pool${pool}/exclude_2309

original_location=$(pwd)
bulk_ref_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/variant_calls/hgsoc-new"
cellsnp_location="$original_location/results/cellSNP/pool${pool}"
output_location="$original_location/results/vireo/${data_type}/pool${pool}/exclude_2309"

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
	bulk_vcf_file=$bulk_ref_location/bulk_diss_polyA_all_samples.vcf.gz
elif [ "${data_type}" == "bulk" ]; then 
	bulk_vcf_file=$bulk_ref_location/bulk_chunk_ribo_all_samples.vcf.gz
else
	echo "Invalid data type specified. Use 'bulk' or 'diss_bulk'."
	exit 1
fi

# Separate the comma-separated string of samples to exclude into an array
IFS=',' read -ra samples <<< "$sample_ids"

# Exclude low quality samples from the bulk VCF file
# Match sample_to_exclude (which is only a partial match) with the sample names in the VCF file and exclude it
for sample in "${samples[@]}"; do
	echo "Excluding sample: $sample"
	full_sample_to_exclude=$(bcftools query -l "$bulk_vcf_file" | grep -E "^HGSOC-(${sample})_" || true)
	full_samples_to_exclude+="$full_sample_to_exclude,"
done
full_samples_to_exclude="${full_samples_to_exclude%,}"  # Remove the trailing comma

# full_sample_to_exclude=$(bcftools query -l "$bulk_vcf_file" | grep -E "^HGSOC-(${samples_to_exclude})_" || true)
# concatenate full samples separate by comma
# full_samples_to_exclude=$(printf "%s\n" "$full_samples_to_exclude" | paste -sd ',' -)

if [ -z "$full_samples_to_exclude" ]; then
	echo "No VCF sample matched samples_to_exclude=${samples_to_exclude}."
	echo "Expected a sample starting with HGSOC-(${samples_to_exclude})_."
	exit 1
fi

match_count=$(printf "%s\n" "$full_sample_to_exclude" | wc -l | tr -d ' ')
if [ "$match_count" -ne 1 ]; then
	echo "Expected exactly 1 VCF sample match for samples_to_exclude=${samples_to_exclude}, found ${match_count}."
	echo "Matches:"
	printf "%s\n" "$full_samples_to_exclude"
	exit 1
fi

filtered_bulk_vcf_file=${bulk_vcf_file%.vcf.gz}_rehead.vcf.gz
bcftools view -s ^${full_samples_to_exclude} "$bulk_vcf_file" -Oz -o "$filtered_bulk_vcf_file" 
bcftools index -t -f "$filtered_bulk_vcf_file"

# Run vireo
vireo \
	-c $cellsnp_location \
	-N $n_samples \
	-o $output_location \
	-d $filtered_bulk_vcf_file \
	--randSeed=12

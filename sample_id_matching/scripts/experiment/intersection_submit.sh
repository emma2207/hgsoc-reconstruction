#!/bin/sh

#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=intersect
#SBATCH --output=intersect_%A_%a.log
#SBATCH --error=intersect_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
module load bcftools
conda activate gnomad-toolbox

sample_ids=(
    "2251"
    "2267"
    "2283"
    "2293"
    "2380"
    "2428"
    "2467"
    "2497"
)

original_location=$(pwd)
data_location="${original_location}/results/1a_individual_vcf/hgsoc/real_data/bulk_dissociated_polyA_vs_single-cell/ncells_null"
output_location="${original_location}/intersection_results"

mkdir -p $output_location

# Find all .vcf.gz files in the data location with the sample ID in their name that
# matches the current SLURM_ARRAY_TASK_ID
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing sample ID: ${sample_id}"

vcf_files=($(find ${data_location} -name "*.vcf.gz" -exec basename {} \; | grep "${sample_id}"))

echo "Found ${#vcf_files[@]} VCF files."

# # Intersect variant data with GnomAD variants of with $sample_id in their name
for file in "${vcf_files[@]}"; do
    python intersection_test.py "${data_location}/${file}"
done

# Index the resulting .vcf.gz files
for file in "${vcf_files[@]}"; do
    zcat ${output_location}/${file} | bgzip -c > ${output_location}/${file}.tmp
    mv ${output_location}/${file}.tmp ${output_location}/${file}
    bcftools index ${output_location}/${file}
done


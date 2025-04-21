#!/bin/sh

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --mem=5G
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools
#SBATCH --output=02_make__diss_bulk_vcf_file_%A_%a.log
#SBATCH --error=02_make_diss_bulk_vcf_file_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load bcftools/1.16

# Based on https://github.com/greenelab/deconvolution_pilot, 
# file scripts/genetic_dmx/02_make_bulk_vcf_file.sh

export pool=$SLURM_ARRAY_TASK_ID
pool_csv=$1

# Script location is /scratch/alpine/$USER/hgsoc/pooled
original_location=`pwd`
index_location="/projects/$USER/hgsoc/refdata-gex-GRCh38-2024-A"
bulk_bam_location="/scratch/alpine/$USER/hgsoc/bulk"
mkdir -p bcftools/pool$pool
output_location="$original_location/bcftools/pool$pool"

# Figure out what samples are in the pool
while IFS=, read -r index col1 col2 col3 col4
do
	if [[ "$index" == "Pool$pool" ]]; then
		col1=${col1%.0}
		col2=${col2%.0}
		col3=${col3%.0}
		col4=${col4%.0}
		sample_array=($col1 $col2 $col3 $col4)
		echo ${sample_array[*]}
	fi
done < $pool_csv

# Check which folder the bulk data for each sample is in
samples_230509=("2018" "2023" "2094" "2126" "2129" \
	"2186" "2202" "2209" "2240" "2309" "2313" "2401" \
	"2407" "2444" "2483" "2526")
samples_230626=("2216" "2221" "2230" "2238" "2246" \
	"2249" "2268" "2278" "2296" "2364" "2408" "2416" \
	"2423" "2430" "2455" "2460" "2466" "2477" "2507" \
	"2514")

bam_files=()
for sample in ${sample_array[@]}
do
	if [[ " ${samples_230509[*]} " =~ " $sample " ]]; then
		bam_files+=("$bulk_bam_location/230509/${sample}/STAR/Aligned.sortedByCoord.out.bam")
	elif [[ " ${samples_230626[*]} " =~ " $sample " ]]; then
		bam_files+=("$bulk_bam_location/230626/${sample}/STAR/Aligned.sortedByCoord.out.bam")
	else
		echo "Sample $sample not found in bulk"
	fi
done

echo ${bam_files[*]}

# Run bcftools
bcftools mpileup -Ou \
	-f $index_location/fasta/genome.fa \
	${bam_files[@]} | \
bcftools call -mv -Ov \
	-o $output_location/bcftools_bulk_pool$pool.vcf

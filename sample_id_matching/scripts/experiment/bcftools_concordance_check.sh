#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools_discordance
#SBATCH --output=bcftools_discordance_%J.log
#SBATCH --error=bcftools_discordance_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load bcftools

dataset="hgsoc"
mod1="bulk_dissociated_polyA"
mod2="single-cell"

# Current location is /scratch/alpine/${USER}/hgsoc/simulations/scripts/experiment
current_location=`pwd`
ref_genome="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
input_location="${current_location}/results/1a_vcf/${dataset}/real_data/${mod1}_vs_${mod2}/ncells_null/read_depth_0"
output_location="${current_location}/test/${dataset}"
mkdir -p "${output_location}"

file1="${input_location}/${mod1}_modality_variants.vcf.gz"
file2="${input_location}/${mod2}_modality_variants.vcf.gz"

file1_norm="${output_location}/file1.norm.vcf.gz"
file2_norm="${output_location}/file2.norm.vcf.gz"
file1_gtcheck="${output_location}/file1.gtcheck.vcf.gz"
file2_gtcheck="${output_location}/file2.gtcheck.vcf.gz"
gtcheck_output="${output_location}/test_output.txt"

bcftools norm -f "${ref_genome}" -m- -Oz -o "${file1_norm}" "${file1}"
bcftools norm -f "${ref_genome}" -m- -Oz -o "${file2_norm}" "${file2}"
bcftools index -f "${file1_norm}"
bcftools index -f "${file2_norm}"

bcftools view -v snps -m2 -M2 -r \
    chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    -Oz -o "${file1_gtcheck}" "${file1_norm}"
bcftools view -v snps -m2 -M2 -r \
    chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    -Oz -o "${file2_gtcheck}" "${file2_norm}"

bcftools index -f "${file1_gtcheck}"
bcftools index -f "${file2_gtcheck}"

bcftools gtcheck -g "${file1_gtcheck}" "${file2_gtcheck}" > "${gtcheck_output}"

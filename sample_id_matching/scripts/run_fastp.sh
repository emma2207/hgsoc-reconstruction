#!/bin/sh

#SBATCH --array=2,5-7
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --ntasks=3
#SBATCH --account=amc-general
#SBATCH --job-name=fastp
#SBATCH --output=fastp_%A_%a.log
#SBATCH --error=fastp_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

# Activate the conda environment for mismatch_project
module load anaconda
conda activate mismatch_project-install

export pool=$SLURM_ARRAY_TASK_ID

sc_fastq_location="/pl/active/cgreene-sc-hgsoc/ariel_sc_HGSOC/GreeneC_OrmanD_scV3_HM5CCDSX5_230201/fastq"
mkdir -p "fastp"
output_location="/scratch/alpine/$USER/hgsoc/pooled/fastp"

# Combine fastq files across lanes
# cat ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L001_R1_001.fastq.gz \
#     ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L002_R1_001.fastq.gz \
#     ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L003_R1_001.fastq.gz \
#     ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L004_R1_001.fastq.gz \
#     > ${output_location}/Pool${pool}-GEX_R1_merged.fastq.gz
# cat ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L001_R2_001.fastq.gz \
#     ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L002_R2_001.fastq.gz \
#     ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L003_R2_001.fastq.gz \
#     ${sc_fastq_location}/Pool${pool}-GEX-CaseyGreene-*_L004_R2_001.fastq.gz \
#     > ${output_location}/Pool${pool}-GEX_R2_merged.fastq.gz

# Run fastp
fastp --in1 ${output_location}/Pool${pool}-GEX_R1_merged.fastq.gz \
      --in2 ${output_location}/Pool${pool}-GEX_R2_merged.fastq.gz \
      --out1 ${output_location}/Pool${pool}_trimmed_R1.fastq.gz \
      --out2 ${output_location}/Pool${pool}_trimmed_R2.fastq.gz \
      --html ${output_location}/Pool${pool}_fastp_report.html \
      --json ${output_location}/Pool${pool}_fastp_report.json \
      --thread 3

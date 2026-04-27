#!/bin/sh

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --ntasks=8
#SBATCH --account=amc-general
#SBATCH --job-name=star
#SBATCH --output=get_cell_by_gene_matrix_%J.log
#SBATCH --error=get_cell_by_gene_matrix_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load star/2.7.10b

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
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID-1]}

# Current location is /scratch/alpine/${USER}/hgsoc/simulations/scripts/experiment
current_location=`pwd`
genome_location="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A"
input_location="/pl/active/cgreene-sc-hgsoc/mismatch_project_data/aligned_reads/hgsoc/single-cell/${sample_id}_single_cell"
output_location="/scratch/alpine/${USER}/hgsoc/simulations/scripts/experiment/miqc/hgsoc"

# Run STARsolo to create cell x gene matrix from Aligned.sortedByCoord.out.bam files
STAR --runThreadN 8 \
     --genomeDir ${genome_location}/star \
     --readFilesType SAM PE \
     --readFilesIn ${input_location}/Aligned.sortedByCoord.out.bam \
     --readFilesCommand "samtools view" \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist None \
     --soloCBlen 16 \
     --soloUMIlen 12 \
     --soloBarcodeReadLength 28 \
     --soloInputSAMattrBarcodeSeq CB UB \
     --soloInputSAMattrBarcodeQual - \
     --soloFeatures Gene \
     --outFileNamePrefix ${output_location}/

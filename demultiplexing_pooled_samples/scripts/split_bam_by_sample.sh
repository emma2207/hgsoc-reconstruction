#!/bin/sh

#SBATCH --array=2-6
#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=split-bam
#SBATCH --output=split-bam_%A_%a.log
#SBATCH --error=split-bam_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load anaconda
conda activate subset-bam-install

export pool=$SLURM_ARRAY_TASK_ID  # 1-10

# Adjust the number of samples for pools 9 & 10
if [ "$pool" -lt 9 ]; then
    nr_samples=4
elif [ "$pool" -lt 10 ]; then
    nr_samples=3
else
    nr_samples=2
fi

# Original location is /scratch/alpine/$USER/hgsoc/pooled
    original_location=`pwd`

# for loop over all samples
for sample in $(seq 1 $nr_samples)
do  
    sample_location="$original_location/bam/pool$pool"
    barcodes_location="$original_location/htseq-count/pool$pool/sample$sample"
    mkdir -p subset-bam/pool$pool/sample$sample

    subset-bam -b $sample_location/pooled.bam \
        -c $barcodes_location/barcodes_per_sample.csv \
        -o "subset-bam/pool$pool/sample$sample/sample.bam" \
        --cores 1
done

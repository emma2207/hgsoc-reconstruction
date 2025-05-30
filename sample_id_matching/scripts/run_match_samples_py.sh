#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=bcftools
#SBATCH --output=match_samples_%J.log
#SBATCH --error=match_samples_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

module load python/3.10.2
module load anaconda
conda activate python_vireosnp_install

# Script location is /scratch/alpine/$USER/hgsoc/pooled

python match_sample_ids.py

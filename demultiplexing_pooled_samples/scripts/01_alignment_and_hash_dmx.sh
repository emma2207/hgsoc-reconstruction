#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --account=amc-general
#SBATCH --job-name=cellranger_multi
#SBATCH --output=01_cellranger_multi_%A_%a.log
#SBATCH --error=01_cellranger_multi_%A_%a.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

set -euo pipefail

module use --append /projects/$USER/lmod-files
module load cellranger/9.0.1
# module load cellranger/7.1.0

# Work around a Cell Ranger 9 telemetry collector crash seen on some HPC runs.
export TENX_DISABLE_TELEMETRY=true

export pool=$SLURM_ARRAY_TASK_ID
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
CSV_PATH="${SUBMIT_DIR}/cellranger_input/cellranger_multi_input_pool${pool}.csv"

if [[ ! -f "${CSV_PATH}" ]]; then
	echo "ERROR: Missing Cell Ranger CSV: ${CSV_PATH}" >&2
	exit 1
fi

LOCALCORES="${SLURM_CPUS_PER_TASK:-16}"
if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" -gt 0 ]]; then
	# SLURM_MEM_PER_NODE is MB; cellranger expects GB for --localmem.
	LOCALMEM_GB=$(( SLURM_MEM_PER_NODE / 1024 - 4 ))
	if [[ "${LOCALMEM_GB}" -lt 8 ]]; then
		LOCALMEM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))
	fi
else
	LOCALMEM_GB=60
fi

# Location is /scratch/alpine/$USER/hgsoc/demultiplexing

cd "${SUBMIT_DIR}"

cellranger multi \
	--id=pool$pool \
	--csv="${CSV_PATH}" \
	--localcores "${LOCALCORES}" \
	--localmem "${LOCALMEM_GB}"

#!/bin/bash

#SBATCH --nodes=1
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --ntasks=20
#SBATCH --account=amc-general
#SBATCH --job-name=cellranger_mkref
#SBATCH --output=00_cellranger_mkref_%J.log
#SBATCH --error=00_cellranger_mkref_%J.err
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

# Build a Cell Ranger 9-compatible reference using cellranger mkref.
#
# IMPORTANT: Do NOT use standalone STAR (e.g., module load star/2.7.10b) to
# rebuild the STAR index for Cell Ranger. STAR 2.7.10b and older write a
# 'genomeType' field in genomeParameters.txt that Cell Ranger 9's bundled
# STAR 2.7.11b does not recognize, causing:
#   FATAL INPUT ERROR: unrecognized parameter name "genomeType"
#
# cellranger mkref uses CR9's own bundled STAR to build the index, which
# guarantees compatibility.
#
# NOTE: 00_create_STAR_index.sh builds a separate standalone STAR index
# used only for bulk RNA-seq alignment (scripts 02+). It is NOT used by
# Cell Ranger and should remain as-is.
#
# After this job completes, update the 'reference' field in all
# cellranger_input/cellranger_multi_input_pool*.csv files to point to
# the new reference path printed at the end of this script.

set -euo pipefail

module purge
module use --append /projects/$USER/lmod-files
module load cellranger/9.0.1

OLD_REF="/projects/${USER}/hgsoc/refdata-gex-GRCh38-2024-A"
NEW_REF_NAME="refdata-gex-GRCh38-2024-A-CR9"
OUTPUT_DIR="/projects/${USER}/hgsoc"

if [[ -d "${OUTPUT_DIR}/${NEW_REF_NAME}" ]]; then
	echo "Reference already exists at ${OUTPUT_DIR}/${NEW_REF_NAME}"
	echo "Delete that directory first if you want to rebuild."
	exit 0
fi

cd "${OUTPUT_DIR}"

cellranger mkref \
	--genome="${NEW_REF_NAME}" \
	--fasta="${OLD_REF}/fasta/genome.fa" \
	--genes="${OLD_REF}/genes/genes.gtf" \
	--nthreads=20

echo ""
echo "Done. New Cell Ranger 9-compatible reference: ${OUTPUT_DIR}/${NEW_REF_NAME}"
echo "Update the 'reference' row in all cellranger_input/cellranger_multi_input_pool*.csv files to:"
echo "  reference,${OUTPUT_DIR}/${NEW_REF_NAME}"

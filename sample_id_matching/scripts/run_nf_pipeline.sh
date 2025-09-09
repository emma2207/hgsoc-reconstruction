#!/bin/bash
#SBATCH --job-name=simulations_nf
#SBATCH --account=amc-general
#SBATCH --output=sim_nf_%J.log
#SBATCH --error=sim_nf_%J.err
#SBATCH --time=01:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --mem=4G
#SBATCH --nodes=1 
#SBATCH --mail-user=emma.lathouwers@cuanschutz.edu
#SBATCH --mail-type=ALL

set -eo pipefail

# --------------------------------------------------
# 0) Resolve project root (directory of this script)
# --------------------------------------------------
# Check if SLURM_SUBMIT_DIR is set (indicates Slurm environment)
if [ -n "${SLURM_SUBMIT_DIR}" ]; then
    # Running on HPC via Slurm
    PRJ_DIR="${SLURM_SUBMIT_DIR}"
    RUN_MODE="HPC"
    echo "Running on HPC (Slurm). Project directory: ${PRJ_DIR}"
else
    # Running locally
    # Get the directory where the script itself is located when run locally
    PRJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    RUN_MODE="LOCAL"
    echo "Running locally. Project directory: ${PRJ_DIR}"
fi

# Change to the project directory
cd "${PRJ_DIR}"
echo "Current working directory: $(pwd)"

# --------------------------------------------------
# 0.5) Specifying paths specific to Alpine
# --------------------------------------------------
# export paths specific to Alpine
if [ "${RUN_MODE}" == "HPC" ]; then
    export PATH=/usr/include:$PATH
    export CPATH=/usr/include/:$CPATH
    export C_INCLUDE_PATH=/usr/include/:$C_INCLUDE_PATH
fi

# --------------------------------------------------
# 1) Load miniforge, activate Conda environment sim_conda
# --------------------------------------------------
# module load miniforge

ENV_YML="${PRJ_DIR}/simulation_environment.yml"

# create the env once; reuse afterwards
 if ! conda env list | grep -q '^simulation_environment '; then
    echo "••• Creating Conda environment simulation_environment from ${ENV_YML}"
    conda env create -f "${ENV_YML}"
 fi

echo "••• Activating Conda environment simulation_environment"
conda init bash
conda activate simulation_environment

# --------------------------------------------------
# 2) Ensure that Nextflow is installed
# --------------------------------------------------
if ! command -v nextflow &> /dev/null
then
    echo "Error: Nextflow is not installed or not found in your PATH."
    echo "Please install Nextflow and ensure it's accessible in your system's PATH."
    echo "You can typically install it with: curl -s https://get.nextflow.io | bash"
    echo "Then move the 'nextflow' executable to a directory in your PATH (e.g., ~/bin or /usr/local/bin)."
    exit 1 # Exit the script with an error code
fi

# --------------------------------------------------
# 3) Run the pipeline
# --------------------------------------------------
echo "••• Launching Nextflow"

# Define Nextflow's work directory based on run mode
NEXTFLOW_WORK_DIR="${PRJ_DIR}/nextflow"

if [ "${RUN_MODE}" == "HPC" ]; then
    nextflow run main.nf -profile slurm -resume -w "${NEXTFLOW_WORK_DIR}" -process.echo
else
    nextflow run main.nf -profile local -resume -w "${NEXTFLOW_WORK_DIR}"
fi

echo "••• Pipeline finished 🎉"

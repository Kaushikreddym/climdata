#!/bin/bash
#SBATCH --job-name=BCSD
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=180G
#SBATCH --time=96:00:00
#SBATCH --output=/beegfs/muduchuru/pkgs_fnl/climdata/usecase/sdba/logs/%x_%j.out
#SBATCH --error=/beegfs/muduchuru/pkgs_fnl/climdata/usecase/sdba/logs/%x_%j.err

#
# SLURM submission script for NEX-GDDP-HYRAS BCSD processing on highmem nodes
# Usage: sbatch submit_highmem.sh <model> <scenario>
# Example: sbatch submit_highmem.sh GFDL-ESM4 ssp370

set -e

# Default values
MODEL=${1:-"GFDL-ESM4"}
SCENARIO=${2:-"ssp370"}

# Validate arguments
if [ -z "$MODEL" ] || [ -z "$SCENARIO" ]; then
    echo "Usage: sbatch submit_highmem.sh <model> <scenario>"
    echo "Example: sbatch submit_highmem.sh GFDL-ESM4 ssp370"
    exit 1
fi


# Create logs directory if it doesn't exist
mkdir -p /beegfs/muduchuru/pkgs_fnl/climdata/usecase/sdba/logs

echo "=========================================="
echo "BCSD Processing Job"
echo "=========================================="
echo "Model    : ${MODEL}"
echo "Scenario : ${SCENARIO}"
echo "Job ID   : ${SLURM_JOB_ID}"
echo "Time     : $(date)"
echo "=========================================="

# Navigate to script directory
cd /beegfs/muduchuru/pkgs_fnl/climdata/usecase/sdba

# Activate conda environment
source /home/muduchuru/.bashrc
conda activate sdba

echo "Python executable: $(which python)"
echo "Python version: $(python --version)"

# Run the Python script with model and scenario as arguments
python -u nexgddp_hyras_bc_zarr.py \
    --model "${MODEL}" \
    --scenarios "${SCENARIO}" \
    --n-processes 40
echo "=========================================="
echo "Job completed at $(date)"
echo "=========================================="

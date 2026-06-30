#!/bin/bash
# Submit all 15 BC jobs: 5 models × 3 scenarios
# Usage: bash submit_all.sh

MODELS=(ACCESS-CM2 CanESM5 EC-Earth3 GFDL-ESM4 MIROC6)
SCENARIOS=(historical ssp370 ssp126)

SCRIPT=/beegfs/muduchuru/pkgs_fnl/climdata/usecase/sdba/submit_highmem.sh
LOGDIR=/beegfs/muduchuru/pkgs_fnl/climdata/usecase/sdba/logs

mkdir -p "${LOGDIR}"

for MODEL in "${MODELS[@]}"; do
    for SCENARIO in "${SCENARIOS[@]}"; do
        JOB_ID=$(sbatch \
            --job-name="BCSD_${MODEL}_${SCENARIO}" \
            "${SCRIPT}" "${MODEL}" "${SCENARIO}" \
            | awk '{print $NF}')
        echo "Submitted ${MODEL} ${SCENARIO}  →  job ${JOB_ID}"
    done
done

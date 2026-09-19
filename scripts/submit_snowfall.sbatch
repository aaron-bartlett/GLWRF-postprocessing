#!/bin/bash --login

#SBATCH -A glrcm
#SBATCH --job-name=postproc_snowfall
#SBATCH --partition=u1-compute
#SBATCH --ntasks=1
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --output=slurm_snowfall_%A_%a.log

# --mem/--time are initial estimates from the data-volume analysis in
# PARALLELIZATION.md, not real profiling. After the first real runs, check
# actual usage with:
#   sacct -j <jobid> --format=JobID,MaxRSS,Elapsed,State
# and adjust the values above.

cd /scratch4/NAGAPE/glrcm/Aaron.Bartlett/POSTPROCESSING
module load rdhpcs-conda
conda activate wrf

source ./scripts/postproc_common.sh

key="${KEYS[$SLURM_ARRAY_TASK_ID]}"
path="${paths[$key]}"
year="${years[$key]}"
scenario="${scenarios[$key]}"

echo "Snowfall Processing for path: $path start_year: $year scenario: $scenario"
python -u snowfall_mean_processing.py \
	-p "$path" -y "$year" -s "$scenario" > "log.snowfall_${key}"

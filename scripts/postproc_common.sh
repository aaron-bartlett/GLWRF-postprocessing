#!/bin/bash
# Shared scenario/path/year lookup tables for the postprocessing SLURM job
# arrays (submit_*.sbatch). Source this file; do not execute it directly.
#
# SLURM_ARRAY_TASK_ID (0-4) indexes KEYS to pick which of the 5 scenarios a
# given array task runs.

declare -a KEYS=("HH" "M2" "M5" "L2" "L5")

declare -A scenarios=(
	["M2"]="ssp245"
	["M5"]="ssp585"
	["L2"]="ssp245"
	["L5"]="ssp585"
	["HH"]="historical"
)
declare -A years=(
	["M2"]="2045"
	["M5"]="2045"
	["L2"]="2085"
	["L5"]="2085"
	["HH"]="2005"
)
declare -A paths=(
	["M2"]="../WRF/run_midcentury_ssp245"
	["M5"]="../WRF/run_midcentury_ssp585"
	["L2"]="../WRF/run_latecentury_ssp245"
	["L5"]="../WRF/run_latecentury_ssp585"
	["HH"]="../WRF/run"
)

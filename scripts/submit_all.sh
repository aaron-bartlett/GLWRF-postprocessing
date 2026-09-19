#!/bin/bash --login
#
# Submits the 6 postprocessing scripts as 6 independent SLURM job arrays
# (submit_snowfall.sbatch, submit_precip.sbatch, submit_monthly.sbatch,
# submit_6hextreme.sbatch, submit_extreme.sbatch, submit_lake.sbatch), each
# with 5 array tasks (one per scenario: HH, M2, M5, L2, L5). That's 30
# independent SLURM jobs total instead of the single ~8-hour serial job this
# script used to submit directly.
#
# Each of the 6 sbatch files has its own --mem/--time sized to that specific
# script (see PARALLELIZATION.md), so a script that needs more memory no
# longer shares a memory ceiling with the other 5, and an OOM/failure in one
# script no longer blocks or kills the others. To resubmit just one script
# (e.g. after tuning its --mem), run its sbatch file directly, e.g.:
#   sbatch submit_snowfall.sbatch


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for sbatch_file in submit_snowfall.sbatch submit_precip.sbatch submit_monthly.sbatch \
		   submit_6hextreme.sbatch submit_extreme.sbatch submit_lake.sbatch; do
	echo "Submitting $sbatch_file (5 array tasks: HH, M2, M5, L2, L5)"
	sbatch "${SCRIPT_DIR}/${sbatch_file}"
done


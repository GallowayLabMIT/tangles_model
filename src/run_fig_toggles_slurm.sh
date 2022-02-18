#!/bin/bash
# Launch with LLsub src/run_fig_toggles_slurm.sh [5,70,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_fig_toggles.jl $LLSUB_RANK $LLSUB_SIZE

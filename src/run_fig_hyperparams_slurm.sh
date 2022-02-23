#!/bin/bash
# Launch with LLsub src/run_fig_hyperparams_slurm.sh [10,60,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_fig_hyperparams.jl $LLSUB_RANK $LLSUB_SIZE

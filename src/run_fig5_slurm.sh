#!/bin/bash
# Launch with LLsub src/run_fig5_slurm.sh [5,70,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_fig5.jl $LLSUB_RANK $LLSUB_SIZE

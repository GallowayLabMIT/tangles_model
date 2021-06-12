#!/bin/bash
# Launch with LLsub src/run_julia_slurm.sh [1,40,1]
source /etc/profile
module load julia/1.5.2

julia --project=. src/simulate_summary.jl $LLSUB_RANK $LLSUB_SIZE

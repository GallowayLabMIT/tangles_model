#!/bin/bash
# Launch with LLsub src/run_simplified_zinani_slurm.sh [5,60,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_simplified_zinani.jl $LLSUB_RANK $LLSUB_SIZE

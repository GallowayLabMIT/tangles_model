#!/bin/bash
# Launch with LLsub src/run_fig_sc_bursting.sh [10,30,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_sc_bursting.jl $LLSUB_RANK $LLSUB_SIZE

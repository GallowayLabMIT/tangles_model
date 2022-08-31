#!/bin/bash
# Launch with LLsub src/run_different_topo_rates_slurm.sh [10,60,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_different_topo_rates.jl $LLSUB_RANK $LLSUB_SIZE

#!/bin/bash
# Launch with LLsub src/run_fig_energy_perturbations_slurm.sh [10,60,1]
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_fig_energy_perturbations.jl $LLSUB_RANK $LLSUB_SIZE

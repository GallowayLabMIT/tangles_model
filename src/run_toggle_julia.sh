#!/bin/bash
# SBATCH -s 1
source /etc/profile
module load julia/1.6.1

julia --project=. src/simulate_toggles.jl

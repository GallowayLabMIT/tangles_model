#!/bin/bash

#SBATCH -o output/stdout/topo_sweep-%j.txt
#SBATCH -c 40
# After initalizing the virtual environment (python -m venv env; pip install -r requirements.txt)
# Run with LLsub run_param_sweep_cluster.sh

# Init modules
source /etc/profile
module load anaconda/2020b
source env/bin/activate

python run_sc_dependent_only.py

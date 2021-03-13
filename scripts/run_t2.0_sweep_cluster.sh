#!/bin/bash

#SBATCH -o output/stdout/t2.0-sweep-%j.txt
#SBATCH -c 40
# After initalizing the virtual environment (python -m venv env; pip install -r requirements.txt)
# Run with LLsub run_param_sweep_cluster.sh

# Init modules
source /etc/profile
module load anaconda/2020b
source env/bin/activate

python t2.0_sweep.py

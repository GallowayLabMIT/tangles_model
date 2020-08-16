#!/bin/bash

#SBATCH -o output/stdout/run-%j-%a.txt
#SBATCH -a 0-9
#SBATCH -c 20
# After initalizing the virtual environment (python -m venv env; pip install -r requirements.txt)
# Run with LLsub run_param_sweep_cluster.sh

# Init modules
source /etc/profile
module load anaconda/2020b
source env/bin/activate

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_COUNT: " $SLURM_ARRAY_TASK_COUNT
python perform_param_sweep.py --runner_id $SLURM_ARRAY_TASK_ID --num_runners $SLURM_ARRAY_TASK_COUNT --n_simulations 5000
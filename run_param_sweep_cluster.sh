#!/bin/sh

#SBATCH -o output/stdout/run-%j-%a.txt
#SBATCH -a 1-4
# After initalizing the virtual environment (python -m venv env; pip install -r requirements.txt)
# Run with LLsub run_param_sweep_cluster.sh

# Init modules
source /etc/profile
module load anaconda3-5.0.1
source env/bin/activate

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
python perform_param_sweep.py $SLURM_ARRAY_TASK_ID
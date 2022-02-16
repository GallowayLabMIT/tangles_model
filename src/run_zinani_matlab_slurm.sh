#!/bin/bash
# Launch with LLsub run_zinani_matlab_slurm.sh [1,10,1]
# from WITHIN the source directory
source /etc/profile

matlab -nodisplay -r 'zinani_unpaired(${LLSUB_RANK}); zinani_paired(${LLSUB_RANK}); exit'

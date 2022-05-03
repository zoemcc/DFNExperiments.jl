#!/bin/bash

# Slurm sbatch options
#SBATCH -o neuralpde_jobarray.log-%A-%a
#SBATCH -a 1-128
#SBATCH -c 1

# Initialize julia path
source /home/gridsan/zmccarthy/.julia_profile

export LOG_DIR="/home/gridsan/zmccarthy/logs/spme_c_e_first"
export SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
export SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT
echo "args:"
echo $SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_TASK_COUNT
echo $LOG_DIR

# Call your script as you would from the command line
julia $DFNEXPERIMENTS_DIR/test/pybamm_generate_spme.jl 

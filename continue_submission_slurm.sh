#!/bin/bash
#
#SBATCH --job-name=qutip_computation_rrobin
#SBATCH --output=res_qutip_computation.txt
#
#SBATCH --ntasks=1
#SBATCH --time=90:00
#SBATCH --mem-per-cpu=100
#
#SBATCH --array=70-139

source activate qutip-env
srun python simu_continue.py $SLURM_ARRAY_TASK_ID

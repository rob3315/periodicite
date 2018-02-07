#!/bin/bash
#
#SBATCH --job-name=qutip_periodic_rrobin
#SBATCH --output=res_qutip_periodic.txt
#
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=500
#
#SBATCH --array=0-19

source activate qutip-env
srun python Periodic_system.py $SLURM_ARRAY_TASK_ID

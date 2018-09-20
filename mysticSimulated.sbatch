#!/bin/bash
#SBATCH --array=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=10
#SBATCH --mem=8G
#SBATCH mem-per-cpu=8G
#SBATCH -o slurmOutMystic/mystic_%A_%a.out
#SBATCH -e slurmOutMystic/mystic_%A_%a.err

module load python

ID=${SLURM_ARRAY_TASK_ID}

srun python mysticSimulated.py ${ID}

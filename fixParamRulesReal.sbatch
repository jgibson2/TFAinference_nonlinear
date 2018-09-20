#!/bin/bash
#SBATCH --array=1-162
#SBATCH --mem=8G
#SBATCH -o slurmOutMystic/mystic_%A_%a.out
#SBATCH -e slurmOutMystic/mystic_%A_%a.err

module load python

ID=${SLURM_ARRAY_TASK_ID}

srun python fixParamRulesReal.py ${ID}

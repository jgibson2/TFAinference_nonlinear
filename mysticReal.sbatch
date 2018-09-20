#!/bin/bash
#SBATCH --array=1
#SBATCH --cpus-per-task=5
#SBATCH --nodes=40
#SBATCH --mem=32G
#SBATCH mem-per-cpu=32G
#SBATCH -o slurmOutMystic/mystic_%A_%a.out
#SBATCH -e slurmOutMystic/mystic_%A_%a.err

module load python

srun python mysticReal.py

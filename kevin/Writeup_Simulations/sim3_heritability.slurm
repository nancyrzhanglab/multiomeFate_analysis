#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=72:00:00

#SBATCH --job-name=sim3
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/multiomeFate/out/kevin/Writeup_Simulations/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/multiomeFate/out/kevin/Writeup_Simulations/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore sim3_heritability.R /home/users/kzlin/kzlinlab/projects/multiomeFate/out/kevin/Writeup_Simulations/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout

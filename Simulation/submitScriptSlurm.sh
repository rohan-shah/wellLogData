#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=5GB
#SBATCH --job-name=wellLogData
trap "echo recieved SIGUSR1;" SIGUSR1;
R CMD BATCH --no-save --no-restore -- script.R script.Rout.$SLURM_JOBID

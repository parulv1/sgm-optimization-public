#!/bin/bash
#### Specify job name
#SBATCH -J sgm
#### Output file
#SBATCH -o <output path>/"%x"_"%j".out
#### Error file
#SBATCH -e <output path>"%x"_"%j".err
#### number of cores 
#SBATCH -n 20
#### Specify queue
#SBATCH --partition=long
#### memory per core
#SBATCH --mem=2G
#### Maximum run time 
#SBATCH --time=2-00:00:00

### You need to change the following paths
export PATH="/home/pverma2/software/miniconda3/bin:$PATH"

source activate spectrome
set -exu

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

nproc --all

which python

python -u ../scripts/sgm_ad_pet.py

[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" -j "$SLURM_JOB_ID"
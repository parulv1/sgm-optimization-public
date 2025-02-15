#!/bin/bash
#### Specify job name
#SBATCH -J sgm_local
#### Output file
#SBATCH -o <output path>/"%x"_"%j".out
#### Error file
#SBATCH -e <error path>/"%x"_"%j".err
#### Number of tasks
#SBATCH --array=1-36
#### number of cores 
#SBATCH -n 2
#### Specify queue
#SBATCH --partition=long
#### memory per core
#SBATCH --mem=2G
#### Maximum run time 
#SBATCH --time=2-00:00:00


export PATH="/home/pverma2/software/miniconda3/bin:$PATH"
source activate spectrome
set -exu


export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


nproc --all

which python 

python -u ../scripts/sgm_fit_local.py $((${SLURM_ARRAY_TASK_ID}-1))

[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" -j "$SLURM_JOB_ID"
#!/bin/bash
#$ -cwd
#### Specify job name
#$ -N sgm_local
#### Output file
#$ -o  <output path>/$JOB_NAME_$TASK_ID.out
#### Error file
#$ -e <error path>/$JOB_NAME_$TASK_ID.err
#### Number of tasks
#$ -t 1-36:1
#### number of cores 
#$ -pe smp 2
#### Specify queue
#$ -q long.q
#### memory per core
#$ -l mem_free=2G
#### Maximum run time 
#$ -l h_rt=336:00:00

#### You need to change the following paths
export PATH="/home/pverma2/software/miniconda3/bin:$PATH"
source activate spectrome
set -exu


export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


nproc --all

which python 

echo "Task id is $SGE_TASK_ID"

python -u ../scripts/sgm_fit_local.py $((${SGE_TASK_ID}-1))

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
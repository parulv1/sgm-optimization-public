#!/bin/bash
#$ -cwd
#### Specify job name
#$ -N sgm_ad_pet_test
#### Output file
#$ -o $JOB_NAME_$JOB_ID.out
#### Error file
#$ -e $JOB_NAME_$JOB_ID.err
#### number of cores 
#$ -pe smp 20
#### Specify queue
#$ -q long.q
#### memory per core
#$ -l mem_free=2G
#### Maximum run time 
#$ -l h_rt=150:00:00

#### You need to change the following paths
export PATH="/wynton/protected/home/rad-wynton-only/pverma2/software/miniconda3/bin:$PATH"
source activate spectrome
set -exu


export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


nproc --all

which python 

python -u ../scripts/sgm_ad_pet.py

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
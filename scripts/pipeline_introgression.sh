#!/bin/bash
#$ -cwd
#$ -V
#$ -t 1-75000
#$ -N ABC
#$ -l h_data=2G,time=20:00:00
#$ -o ./logs/
#$ -e ./logs/

# setup
. /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
set -e
set -u

# ensure the output directory exists
outdir=$scratch/fwdsims/introgression
mkdir -p ${outdir}

NUM_SIMS_PER_JOB=10

for i in `seq 1 ${NUM_SIMS_PER_JOB}`;
do
  # 1. draw parameters
  m=`python -c "import random; print(random.uniform(0.0,0.1))"`
  talpha=`python -c "import random; print(random.uniform(1,10))"` #will need to change depending on scaling
  # 2. simulate and 3. compute stats
  # note this is untested -- we need to save the stats and the parameters together
  slim -d m=${m} -d talpha=${talpha} scripts/introgression.slim | scripts/compute_stats.py | paste - ${m} ${talpha} >> ${outdir}/${SGE_JOB_ID}.${SGE_TASK_ID}.txt
done

#!/bin/bash
#SBATCH --account=windse
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --job-name=many_small_jobs

 
# Load and activate your conda environment
module purge
module load conda
source activate windse
 
# This stops the solver from filling all 36 cores
# Not using this may significantly slow down execution (100 times slower on Eagle)
export OMP_NUM_THREADS=1


# Number of times all processors run a job = # of jobs/ # of cores per node
N_LOOP=2
# Run 1 job per task
N_JOB_PER_LOOP=$SLURM_NTASKS

# Where to store Logs
mkdir -p LOGS


# Run many small jobs
for((j=0;j<$N_LOOP;j++))
do

  for((i=0;i<$N_JOB_PER_LOOP;i++))
  do

    jobid=$((i + j*N_JOB_PER_LOOP))
    windse run test.yaml -p 'general':'name':test_$jobid > LOGS/log$jobid &

  done

  #Wait until all jobs running on the node are done before launching the next batch
  wait

done

# Finalize
echo "All done"

 

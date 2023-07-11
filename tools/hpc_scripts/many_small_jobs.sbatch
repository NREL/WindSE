#!/bin/bash
#SBATCH --account=windse
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --job-name=many_small_jobs

# This SBATCH file demonstrates how a single parameter file, "test.yaml",
# can be used as a baseline for iteratively testing different simulation options.
# Different options can be passed through the command-line specifiers,
# e.g., -p:'general':'name':$jobname, and run in parallel using as many
# simultaneous instances as the HPC node can manage.  This example file further
# demonstrates how multiple batches of such parallelism can be processed using the
# double FOR loop construction. To submit this job to the scheduler, simply call
# 
# sbatch many_small_jobs.sbatch 
#
# from your HPC connection, replacing the account name, job time, node number, 
# conda environment name, etc., in this file as needed.
 
# Load and activate your conda environment
module purge
module load conda
source activate windse
 
# This stops the solver from filling all 36 cores of the node
# Not using this may significantly slow down execution (100 times slower on Eagle)
export OMP_NUM_THREADS=1


# Run 1 job per core
N_JOB_PER_BATCH=$SLURM_NTASKS
# Number of times all processors run a job = # of jobs/ # job per batch
N_BATCH=2

# Where to store logs
mkdir -p LOGS


# Run many small jobs
for((j=0;j<$N_BATCH;j++))
do

  for((i=0;i<$N_JOB_PER_BATCH;i++))
  do
    # Construct a unique job ID
    jobid=$((i + j*N_JOB_PER_BATCH))

    # Create a unique name for this job and its output
    jobname=$(printf "test_job_%03d" $jobid)
    logname=$(printf "log_%03d" $jobid)

    # Run each as a 1-core job, where -n XX can be adjusted as needed
    # (e.g., srun -n 2... where N_JOB_PER_BATCH=18)
    srun -n 1 windse run test.yaml -p 'general':'name':$jobname > LOGS/$logname &

  done

  #Wait until all jobs running on the node are done before launching the next batch
  wait

done

# Finalize
echo "All done"

 

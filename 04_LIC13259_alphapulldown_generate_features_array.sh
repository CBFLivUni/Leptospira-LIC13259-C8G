#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J alphapulldown
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o gen_features_alphapulldown.%u.%N.%j.out
# Define a standard error file
#SBATCH -e gen_features_alphapulldown.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 1
# Request the number of cores
#SBATCH -n 4
# Specify time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 0-6:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores 
##SBATCH --mem-per-cpu=9000M
# Define job array
#SBATCH --array=1-20
#
# Set the maximum stack size to unlimited
ulimit -s unlimited
# Set OpenMP thread number
export OMP_NUM_THREADS=$SLURM_NTASKS

# Load tensorflow and relevant modules
module purge
module load apps/anaconda3/5.2.0-python37
#use source activate gpu to get the gpu virtual environment
source activate AlphaPulldown
# List all modules
module list
#
#
echo =========================================================   
echo SLURM job: submitted  date = `date`      
date_start=`date +%s`

echo -------------  
echo Job output begins                                           
echo -----------------                                           
echo

hostname

echo "Print the following environmetal variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job user                     : $SLURM_JOB_USER"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"
echo "Submit directory             : $SLURM_SUBMIT_DIR"
echo "Temporary directory          : $TMPDIR"
echo "Submit host                  : $SLURM_SUBMIT_HOST"
echo "Queue/Partition name         : $SLURM_JOB_PARTITION"
echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Hostname of 1st node         : $HOSTNAME"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of processes          : $SLURM_NTASKS"
echo "Number of processes per node : $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"

echo   
echo "Running job:"
echo   

echo =========================================================   

# Script to generate features for baits and candidates for downstream in alphapulldown

# The length of the array job is the number of baits x candidates
# This could be be passed to this script on the command line with the following commands:
# baits=`grep -c "" baits.txt`
# candidates=`grep -c "" candidates.txt`
# count=$(( $baits + $candidates ))

# Make output folder for alphapulldown
mkdir -p alphapulldown/LIC13259_no_his

# Copy the msa files to the alphapulldown folder
cp pdbs/LIC13259_no_his/msa/*.a3m alphapulldown/LIC13259_no_his
cp pdbs/C8/msa/*.a3m alphapulldown/LIC13259_no_his

create_individual_features.py \
  --fasta_paths=LIC13259_no_his_C8_fasta.fasta \
  --data_dir=/mnt/lustre/users/ejohn10/alphafold_database \
  --save_msa_files=True \
  --output_dir=alphapulldown/LIC13259_no_his \
  --max_template_date=2050-01-01 \
  --skip_existing=True \
  --use_mmseqs2=True \
  --db_preset=reduced_dbs \
  --seq_index=$SLURM_ARRAY_TASK_ID



echo =========================================================   
# the ret flag is the return code, so you can spot easily if your code failed.
ret=$?

echo   
echo ---------------                                           
echo Job output ends                                           
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================   
echo SLURM job: finished   date = `date`   
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================   
exit $ret
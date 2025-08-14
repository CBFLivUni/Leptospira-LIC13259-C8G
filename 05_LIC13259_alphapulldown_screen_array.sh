#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job
#SBATCH --export=ALL
# Define job name
#SBATCH -J colab_GPU
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
# %A will be replaced by the value of SLURM_ARRAY_JOB_ID and %a will be replaced by the value of SLURM_ARRAY_TASK_ID.
#SBATCH -o alphapulldown.%u.%N.%j.%A_%a.out
# Define a standard error file
#SBATCH -e alphapulldown.%u.%N.%j.%A_%a.err
# Request the GPU partition (gpu). We don't recommend requesting multiple partitions, as the specifications of the nodes in these partitions are different.
#SBATCH -p gpu
# Request the number of nodes
#SBATCH -N 1
# Request the number of GPUs per node to be used (if more than 1 GPU per node is required, change 1 into Ngpu, where Ngpu=2,3,4)
#SBATCH --gres=gpu:1
# Request the number of CPU cores. (There are 24 CPU cores and 4 GPUs on each GPU node in partition gpu,
# so please request 6*Ngpu CPU cores, i.e., 6 CPU cores for 1 GPU, 12 CPU cores for 2 GPUs, and so on.)
#SBATCH -n 6
# Set time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores
##SBATCH --mem=90GB or #SBATCH --mem-per-cpu=9000M
# Define job array
#SBATCH --array=1-84
#
# Set the maximum stack size to unlimited
ulimit -s unlimited

# Load relevant modules
module purge
#module load apps/anaconda3/5.2.0-python37
micromamba activate AlphaPulldown
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/data1/users/ejohn10/y/envs/AlphaPulldown/

#module unload libs/nvidia-cuda/8.0.61/bin
#module load libs/cuda/11.3/bin
#module unload libs/cudnn/7.0.5_cuda8
#module load libs/cudnn/8.2.1_cuda11.3


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

# Script to carry out alphapulldown analysis between baits and candidates

# The length of the array job is the number of baits x candidates
# This could be be passed to this script on the command line with the following commands:
# baits=`grep -c "" baits.txt`
# candidates=`grep -c "" candidates.txt`
# count=$(( $baits * $candidates ))

run_multimer_jobs.py \
  --mode=pulldown \
  --models_to_relax=all \
  --num_cycle=3 \
  --num_predictions_per_model=10 \
  --output_path=alphapulldown/LIC13259_no_his \
  --data_dir=/mnt/lustre/users/ejohn10/alphafold_database \
  --protein_lists=LIC13259_no_his_baits.txt,C8_all_species_candidates.txt \
  --monomer_objects_dir=alphapulldown/LIC13259_no_his \
  --compress_result_pickles=True \
  --remove_result_pickles=True \
  --job_index=$SLURM_ARRAY_TASK_ID


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
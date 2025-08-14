#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J process_files
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o process_files.%u.%N.%j.out
# Define a standard error file
#SBATCH -e process_files.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 1
# Request the number of cores
#SBATCH -n 6
# Specify time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores 
##SBATCH --mem-per-cpu=9000M
#
# Set the maximum stack size to unlimited
ulimit -s unlimited
# Set OpenMP thread number
export OMP_NUM_THREADS=$SLURM_NTASKS

# Load your own modules
module purge

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
echo "Submit host                  : $SLURM_SUBMIT_HOST"./processed_fasta/ompl1/
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


# Script to carry out multimer clustering
# Should be ran from the project parent directory

# LIC13259 
# Create directories
mkdir -p results/LIC13259_no_his/multimer_clustering

for dir in ./results/LIC13259_no_his/alphapulldown/Leptospira*/ ; do
    if [[ -d "$dir" ]]; then
        NAME=$(basename "${dir}")
        echo "Clustering for $NAME"
        
        # Set up output paths
        OUTPUT_CLUSTER_NAME="${NAME}_cluster"
        OUTPUT_CLUSTER_DIR="./results/LIC13259_no_his/multimer_clustering/"
        
        # Run Foldseek multimer clustering
        foldseek easy-multimercluster "$dir" "$OUTPUT_CLUSTER_DIR/$OUTPUT_CLUSTER_NAME" tmp \
            --multimer-tm-threshold 0.5 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.5
        
        # Move the resulting .tsv file to multimer_clustering directory
        if [[ -f "${OUTPUT_CLUSTER_DIR}/${OUTPUT_CLUSTER_NAME}_cluster.tsv" ]]; then
            mv "${OUTPUT_CLUSTER_DIR}/${OUTPUT_CLUSTER_NAME}_cluster.tsv" "${OUTPUT_CLUSTER_DIR}/${NAME}_cluster_results.tsv"
        else
            echo "Warning: Expected output .tsv file not found for $NAME"
        fi
    fi
done




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
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
# Insert your own username to get e-mail notifications (note: keep just one "#" before SBATCH)
#SBATCH --mail-user=ejohn10@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
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
# Should be run from the project parent directory

# Define source and target directories
SOURCE_DIR="./results/LIC13259_no_his/alphapulldown/C8G"
TARGET_DIR="./results/LIC13259_no_his/alphapulldown/C8G_representative_clusters"
FILE_LIST="LIC13259_no_his_representive_clusters.txt"

# Create the target directory if it doesn't already exist
mkdir -p "$TARGET_DIR"

# Use xargs to read each line from the file list and copy matching files from SOURCE_DIR to TARGET_DIR
xargs -a "$FILE_LIST" -I {} cp "${SOURCE_DIR}/{}" "$TARGET_DIR"

# Define output directory for Foldseek results
OUTPUT_CLUSTER_DIR="./results/LIC13259_no_his/multimer_clustering/selected_clusters_pairwise_comparisons"
mkdir -p "${OUTPUT_CLUSTER_DIR}"

# Run Foldseek multimer easy search on each .pdb file in TARGET_DIR
for pdb in "$TARGET_DIR"/*.pdb; do
  # Check if the .pdb file exists (important for empty directories)
  if [[ -f "$pdb" ]]; then
    NAME=$(basename "$pdb" .pdb)
    QUERY_PDB="$pdb"
    
    echo "Querying ${NAME}.pdb for structural similarity"
    
    # Run Foldseek for multimer clustering
    foldseek easy-multimersearch \
      "$QUERY_PDB" \
      "$TARGET_DIR" "${OUTPUT_CLUSTER_DIR}/${NAME}_structural_similarity.m8" tmp \
      --format-output "query,target,alntmscore,lddt,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits" \
      --exhaustive-search
  else
    echo "No .pdb files found in ${TARGET_DIR} to process."
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
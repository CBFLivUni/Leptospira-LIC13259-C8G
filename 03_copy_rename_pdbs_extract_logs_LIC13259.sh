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


# Script to extract and rename PDB files 
# It also extracts the relevant log information and sequencing depth image
# It only extracts the PDB for the best ranked structure 
# Should be ran from the project parent directory

# Create directories
mkdir -p pdbs/LIC13259_no_his/sequencing_depth/
mkdir -p pdbs/LIC13259_no_his/logs/
mkdir -p pdbs/LIC13259_no_his/jsons/
mkdir -p pdbs/LIC13259_no_his/msa/

# LIC13259 no his
for i in ./processed_fasta/LIC13259_no_his/*/ ; do
    
    NAME=$(basename ${i} /)
    echo "$NAME"
    cp ./processed_fasta/LIC13259_no_his/$NAME/*_relaxed_rank_001_*.pdb ./pdbs/LIC13259_no_his/$NAME.pdb
    cp ./processed_fasta/LIC13259_no_his/$NAME/*coverage.png ./pdbs/LIC13259_no_his/sequencing_depth/"$NAME"_coverage.png
    cp ./processed_fasta/LIC13259_no_his/$NAME/*_rank_001_*.json ./pdbs/LIC13259_no_his/jsons/"$NAME".json
    cp ./processed_fasta/LIC13259_no_his/$NAME/*.a3m ./pdbs/LIC13259_no_his/msa/"$NAME".a3m
    grep 'rank_001' ./processed_fasta/LIC13259_no_his/$NAME/log.txt  > ./pdbs/LIC13259_no_his/logs/TMP.txt
    echo -n "sequencing depth=" $(tr -cd '>' < ./processed_fasta/LIC13259_no_his/$NAME/*.a3m | wc -c)  >> ./pdbs/LIC13259_no_his/logs/TMP.txt
    paste -sd ' ' ./pdbs/LIC13259_no_his/logs/TMP.txt > ./pdbs/LIC13259_no_his/logs/$NAME.txt
    rm ./pdbs/LIC13259_no_his/logs/TMP.txt
    
done


# C8
for i in ./processed_fasta/C8/*/ ; do
    
    NAME=$(basename ${i} /)
    echo "$NAME"
    cp ./processed_fasta/C8/$NAME/*_relaxed_rank_001_*.pdb ./pdbs/C8/$NAME.pdb
    cp ./processed_fasta/C8/$NAME/*coverage.png ./pdbs/C8/sequencing_depth/"$NAME"_coverage.png
    cp ./processed_fasta/C8/$NAME/*_rank_001_*.json ./pdbs/C8/jsons/"$NAME".json
    cp ./processed_fasta/C8/$NAME/*.a3m ./pdbs/C8/msa/"$NAME".a3m
    grep 'rank_001' ./processed_fasta/C8/$NAME/log.txt  > ./pdbs/C8/logs/TMP.txt
    echo -n "sequencing depth=" $(tr -cd '>' < ./processed_fasta/C8/$NAME/*.a3m | wc -c)  >> ./pdbs/C8/logs/TMP.txt
    paste -sd ' ' ./pdbs/C8/logs/TMP.txt > ./pdbs/C8/logs/$NAME.txt
    rm ./pdbs/C8/logs/TMP.txt
    
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
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


# Script to extract preliminary alphapulldown results for ompl1

# It firstly extracts the top ranked .pdb file for that binding, relevant log information and sequencing depth
# Next it extracts all .pdb files for each binding pair combination
# Then it removes .pdb files that have more than one backbone clash using the check_clashes.py script
# Finally, it examines the interfaces identifying:
# -- amino acids most common in the LIC13259 orthologue involved in binding
# -- amino acid binding pairs between LIC13259 and the relevant C8 subunit

# Should be ran from the project parent directory

# Create directories
mkdir -p results/LIC13259_no_his/alphapulldown/model_rankings
mkdir -p results/LIC13259_no_his/alphapulldown/PAE
mkdir -p results/LIC13259_no_his/alphapulldown/contacts_summary
mkdir -p results/LIC13259_no_his/alphapulldown/residue_pairs


# Copy top-ranked .pdb file, log information, and PAE plot
for dir in ./alphapulldown/LIC13259_no_his/*/ ; do
    NAME=$(basename "${dir}")
    echo "$NAME"
    cp "${dir}ranked_0.pdb" "results/LIC13259_no_his/alphapulldown/${NAME}.pdb"
    cp "${dir}"*PAE_plot_ranked_0.png "results/LIC13259_no_his/alphapulldown/PAE/${NAME}_coverage.png"
    sed -n 3,52p "${dir}"ranking_debug.json > "results/LIC13259_no_his/alphapulldown/model_rankings/${NAME}.txt"
done


# Copy all .pdb files for each binding pair combination
for dir in ./alphapulldown/LIC13259_no_his/*/ ; do
    NAME=$(basename "${dir}")
    echo "Processing directory: $NAME"
    DEST_DIR="./results/LIC13259_no_his/alphapulldown/$NAME"
    mkdir -p "$DEST_DIR"

    for file in "${dir}"ranked*.pdb ; do
        FILENAME=$(basename "${file}")
        echo "Copying file: $FILENAME"
        cp -i "${file}" "$DEST_DIR/${FILENAME%.*}_$NAME.pdb"
    done
done


# Identify and remove .pdb files with more than one backbone clash
base_dir="./results/LIC13259_no_his/alphapulldown/"

for dir in "$base_dir"*/ ; do
    echo "Running check_clashes.py script on $dir directory"
    python check_clashes.py "$dir" --output_csv "$dir/clashes_summary.csv" --threshold 3.6

    if [[ -f "$dir/clashes_summary.csv" ]]; then
        while IFS=, read -r pdb_file clashes_count clashes_per_100_residues; do
            [[ "$pdb_file" == "PDB_File" ]] && continue
            if (( clashes_count > 0 )); then
                pdb_path="$dir/$pdb_file"
                if [[ -f "$pdb_path" ]]; then
                    echo "Deleting $pdb_path with $clashes_count clashes"
                    rm "$pdb_path"
                fi
            fi
        done < "$dir/clashes_summary.csv"
    else
        echo "No clashes_summary.csv file found in $dir"
    fi
done


# Rank the .pdb files, so only the top 25 are remaining (after removing clashes)
for dir in "$base_dir"/*; do
    if [[ -d "$dir" ]]; then
        echo "Processing directory: $dir"
        
        # Find all .pdb files with the prefix 'ranked_' in the subdirectory
        pdb_files=("$dir"/ranked_*.pdb)
        
        # Sort the files by their rank number extracted from the filename
        sorted_files=($(for file in "${pdb_files[@]}"; do
            rank=$(echo "$file" | grep -oP 'ranked_\K[0-9]+')
            echo "$rank $file"
        done | sort -n | cut -d' ' -f2-))
        
        # Keep only the top 25 files, delete the rest
        for i in "${!sorted_files[@]}"; do
            if [[ $i -ge 25 ]]; then
                echo "Deleting ${sorted_files[$i]}"
                rm "${sorted_files[$i]}"
            fi
        done
    fi
done


# Identify conserved residues at the interface between LIC13259_no_his and plasminogen for remaining .pdb files
for dir in "$base_dir"*/ ; do
    if [[ -d "$dir" ]]; then
        python find_interface_residues.py "$dir"
    fi
done



# Summarised results for LIC13259_no_his chain
for dir in ./results/LIC13259_no_his/alphapulldown/*/ ; do
    NAME=$(basename "${dir}")
    echo "$NAME"
    cp "${dir}contacts_csv/residue1_summary.csv" "results/LIC13259_no_his/alphapulldown/contacts_summary/${NAME}_LIC13259_no_his_contacts_summary.csv"
done


# Summarised results for specific contacts between LIC13259_no_his and cfh
for dir in ./results/LIC13259_no_his/alphapulldown/*/ ; do
    NAME=$(basename "${dir}")
    echo "$NAME"
    cp "${dir}contacts_csv/summary_contacts.csv" "results/LIC13259_no_his/alphapulldown/residue_pairs/${NAME}_LIC13259_no_his_cfh_residue_pairs.csv"
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
#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J af_db
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o af_db.%u.%N.%j.out
# Define a standard error file
#SBATCH -e af_db.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 1
# Request the number of cores
#SBATCH -n 4
# Specify time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores 
##SBATCH --mem-per-cpu=9000M
# Insert your own username to get e-mail notifications (note: keep just one "#" before SBATCH)
#SBATCH --mail-user=ejohn10@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL

#------------
# Params
# https://github.com/google-deepmind/alphafold/blob/main/scripts/download_alphafold_params.sh


DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/params"
SOURCE_URL="https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
#aria2c "${SOURCE_URL}" --dir="${ROOT_DIR}"

tar --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
  --directory="${ROOT_DIR}" --preserve-permissions
rm "${ROOT_DIR}/${BASENAME}"


https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar

#------------
# PDB mmcif
# https://github.com/google-deepmind/alphafold/blob/main/scripts/download_pdb_mmcif.sh


DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/pdb_mmcif"
RAW_DIR="${ROOT_DIR}/raw"
MMCIF_DIR="${ROOT_DIR}/mmcif_files"

echo "Running rsync to fetch all mmCIF files (note that the rsync progress estimate might be inaccurate)..."
echo "If the download speed is too slow, try changing the mirror to:"
echo "  * rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ (Europe)"
echo "  * ftp.pdbj.org::ftp_data/structures/divided/mmCIF/ (Asia)"
echo "or see https://www.wwpdb.org/ftp/pdb-ftp-sites for more download options."
mkdir --parents "${RAW_DIR}"
rsync --recursive --links --perms --times --compress --info=progress2 --delete --port=33444 \
  rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ \
  "${RAW_DIR}"

echo "Unzipping all mmCIF files..."
find "${RAW_DIR}/" -type f -iname "*.gz" -exec gunzip {} +

echo "Flattening all mmCIF files..."
mkdir --parents "${MMCIF_DIR}"
find "${RAW_DIR}" -type d -empty -delete  # Delete empty directories.
for subdir in "${RAW_DIR}"/*; do
  mv "${subdir}/"*.cif "${MMCIF_DIR}"
done

# Delete empty download directory structure.
find "${RAW_DIR}" -type d -empty -delete

#aria2c "https://files.wwpdb.org/pub/pdb/data/status/obsolete.dat" --dir="${ROOT_DIR}"
wget -O "${ROOT_DIR}/obsolete.dat" "https://files.wwpdb.org/pub/pdb/data/status/obsolete.dat"

#--------------------------
# PDB70
# https://github.com/google-deepmind/alphafold/blob/main/scripts/download_pdb70.sh

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/pdb70"
SOURCE_URL="http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/pdb70_from_mmcif_200401.tar.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
tar --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
  --directory="${ROOT_DIR}"
rm "${ROOT_DIR}/${BASENAME}"


#------------
# uniref90

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/uniref90"
SOURCE_URL="https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${BASENAME}"
popd


#------------
#small_bfd

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/small_bfd"
SOURCE_URL="https://storage.googleapis.com/alphafold-databases/reduced_dbs/bfd-first_non_consensus_sequences.fasta.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${BASENAME}"
popd

#------------
# bfd
DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/bfd"
# Mirror of:
# https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz.
SOURCE_URL="https://storage.googleapis.com/alphafold-databases/casp14_versions/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
tar --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
  --directory="${ROOT_DIR}"
rm "${ROOT_DIR}/${BASENAME}"



#------------
# mgnify
DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/mgnify"
# Mirror of:
# https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_05/mgy_clusters.fa.gz
SOURCE_URL="https://storage.googleapis.com/alphafold-databases/v2.3/mgy_clusters_2022_05.fa.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${BASENAME}"
# gzip -d mgy_clusters_2022_05.fa.gz
popd

#------------
# uniprot
DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/uniprot"

TREMBL_SOURCE_URL="https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
TREMBL_BASENAME=$(basename "${TREMBL_SOURCE_URL}")
TREMBL_UNZIPPED_BASENAME="${TREMBL_BASENAME%.gz}"

SPROT_SOURCE_URL="https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
SPROT_BASENAME=$(basename "${SPROT_SOURCE_URL}")
SPROT_UNZIPPED_BASENAME="${SPROT_BASENAME%.gz}"

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${TREMBL_BASENAME}" "${TREMBL_SOURCE_URL}"
wget -O "${ROOT_DIR}/${SPROT_BASENAME}" "${SPROT_SOURCE_URL}"
pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${TREMBL_BASENAME}"
gunzip "${ROOT_DIR}/${SPROT_BASENAME}"

# Concatenate TrEMBL and SwissProt, rename to uniprot and clean up.
cat "${ROOT_DIR}/${SPROT_UNZIPPED_BASENAME}" >> "${ROOT_DIR}/${TREMBL_UNZIPPED_BASENAME}"
mv "${ROOT_DIR}/${TREMBL_UNZIPPED_BASENAME}" "${ROOT_DIR}/uniprot.fasta"
rm "${ROOT_DIR}/${SPROT_UNZIPPED_BASENAME}"
popd

#gzip -d uniprot_sprot.fasta.gz
#gzip -d uniprot_trembl.fasta.gz


#------------
# Uniref30

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/uniref30"
# Mirror of:
# https://wwwuser.gwdg.de/~compbiol/uniclust/2021_03/UniRef30_2021_03.tar.gz
SOURCE_URL="https://storage.googleapis.com/alphafold-databases/v2.3/UniRef30_2021_03.tar.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
tar --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
  --directory="${ROOT_DIR}"
rm "${ROOT_DIR}/${BASENAME}"


#------------
# Uniclust30

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/uniclust30"
# Mirror of:
# https://wwwuser.gwdg.de/~compbiol/uniclust/2021_03/UniRef30_2021_03.tar.gz
SOURCE_URL="https://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
wget -O "${ROOT_DIR}/${BASENAME}" "${SOURCE_URL}"
tar --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
  --directory="${ROOT_DIR}"
rm "${ROOT_DIR}/${BASENAME}"



date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================
echo SLURM job: finished   date = `date`
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo ==========
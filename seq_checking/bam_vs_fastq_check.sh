#!/bin/bash

set -e
set -o pipefail

#   Define usage message
function usage() {
    echo -e "\
$0: \n\
\n\
This script calls on a Python script to verify the correctness of SAM/BAM file contents.
It does this by comparing Illumina sequence identifiers present in the SAM/BAM file to
sequence identifiers in the FASTQ files. \n\
\n\
Usage: ./bam_vs_fastq_check [SAM/BAM_file_list] [accession_names_list] [fastq_file_list] [script_dir] [out_dir] [scratch_dir]\n\
\n\
NOTE: arguments must be provided in this order! \n\
\n\
where: \n\
1. [SAM/BAM_file_list]  \n\
2. [accession_names_list] \n\
3. [fastq_file_list] \n\
4. [script_dir] \n\
5. [out_dir] \n\
5. [scratch_dir] \n\

Dependencies:
- Python3: >v3.7
- Samtools: >v1.8.0
" >&2
exit 1
}

if [[ $# -lt 1 ]]; then usage; fi

# User provided input arguments
# List containing full filepaths to SAM or BAM files
FILE_LIST=$1
# List of accession names only. Must match part of SAM/BAM and FASTQ filename.
ACC_LIST=$2
# List of FASTQ files that correspond to SAM/BAM samples
FASTQ_LIST=$3
# Full filepath to directory that contains our scripts
# Note: Please make sure bam_vs_fastq_check.sh and check_seq_id.py are in the same directory.
SCRIPT_DIR=$4
# Full filepath to directory to output files
OUT_DIR=$5
# Full filepath to scratch directory to temporarily store intermediate files
#   Note: This can be the same as the out directory but can be useful if you are running short
#   on storage space but have a scratch directory that doesn't count toward your storage.
SCRATCH_DIR=$6

# Check if out directories exists, if not create them
mkdir -p "${OUT_DIR}"
mkdir -p "${SCRATCH_DIR}"
mkdir -p "${SCRATCH_DIR}"/intermediates
# Setup data structures
ACC_ARRAY=($(cat "${ACC_LIST}"))

function compare_seq_id() {
    local accession=$1
    local aligned_list=$2
    local fastq_list=$3
    local out_dir=$4
    # Pull out aligned sample we are currently working with
    aligned=$(awk -v pat="${accession}" '$1 ~ pat {print}' ${aligned_list})
    # Pull out filepaths for forward and reverse reads for sample we are currently working with
    fastq_fwd=$(awk -v pat="${accession}" '$1 ~ pat {print}' ${fastq_list} | grep "R1")
    fastq_rev=$(awk -v pat="${accession}" '$1 ~ pat {print}' ${fastq_list} | grep "R2")

    # Pull out header lines from aligned file and store in file
    # Note: We are not working with the entire SAM/BAM file but instead pulling out the header
    # lines and Illumina sequence identifiers to significantly reduce the amount of data we
    # have to work with.
    samtools view -H "${aligned}" > "${SCRATCH_DIR}"/intermediates/"${accession}"_seqIDs.txt
    # Pull out Illumina sequence identifiers from aligned file and store in a file
    samtools view "${aligned}" | awk '{print $1}' > "${SCRATCH_DIR}"/intermediates/"${accession}"_seqIDs.txt

    # Call on Python script to make comparisons
    # [Placeholder]
    python3 check_seq_id.py
}

export -f compare_seq_id

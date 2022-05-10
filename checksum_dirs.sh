#!/bin/bash

set -e
set -o pipefail

# This script generates checksums for files in two directories on MSI
# Use case: compare original files to transferred files

# Dependencies
module load parallel/20210822

# User provided command line arguments
DIR1="$1"
DIR2="$2"
# Directory we want to output the checksum files to
OUT_DIR="$3"
FILE_SUFFIX="$4"

function usage() {
    echo -e "\
Usage: `basename $0` [dir1] [dir2] [out_dir] [file_suffix]
Where:
    1) [dir1] is the full filepath to directory 1 of files
    2) [dir2] is the full filepath to directory 2 of files
    3) [out_dir] is the directory we want to store the checksum files
    4) [file_suffix] is the file suffix that matches files in dir1/dir2 we
            want to generate checksums for.
" >&2
    exit 1
}

export -f usage

# Check if we have less than one argument
#   If yes, display usage message
if [[ "$#" -lt 1 ]]; then usage; fi

# Checksum on multiple files in parallel
# Run checksums on DIR1
echo "Generating checksums for directory 1..."
cd ${DIR1}
find *"${FILE_SUFFIX}" | parallel "md5sum {}" > ${OUT_DIR}/md5_dir1_unsorted.txt
# Run checksums on DIR2
echo "Generating checksums for directory 2..."
cd ${DIR2}
find *"${FILE_SUFFIX}" | parallel "md5sum {}" > ${OUT_DIR}/md5_dir2_unsorted.txt

# Sort checksum files
echo "Sorting checksum files..."
sort -k 2,2 ${OUT_DIR}/md5_dir1_unsorted.txt > ${OUT_DIR}/md5_dir1.txt
sort -k 2,2 ${OUT_DIR}/md5_dir2_unsorted.txt > ${OUT_DIR}/md5_dir2.txt

# Cleanup intermediate files
rm ${OUT_DIR}/md5_dir1_unsorted.txt
rm ${OUT_DIR}/md5_dir2_unsorted.txt

echo "Done."

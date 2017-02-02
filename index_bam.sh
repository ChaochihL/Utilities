#!/bin/bash

#PBS -l mem=1gb,nodes=1:ppn=8,walltime=6:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -u
set -o pipefail

#   Define usage message
function usage() {
    echo -e "\
$0: \n\
\n\
Usage: ./index_bam.sh [bamFile_list.txt] [bamFile_dir] \n\
\n\
where: \n\
1. [bamFile_list] is a list of bam files to be indexed
2. [bamFile_dir] is where our BAM files are located
" >&2
exit 1
}

if [[ $# -lt 3 ]]; then usage; fi

#   Dependencies
module load samtools
module load parallel

#   Function to index bam files
function indexBAM() {
    local BAMFile="$1"
    local out_dir="$2"
    #   Sample name
    sampleName=`basename "${BAMFile}" .bam`
    #   Create BAI index for BAM file
    samtools index -b "${out_dir}/Finished/${sampleName}.bam"
}

#   Export function
export -f indexBAM

#   Arguments provided by user
BAM_LIST=$1 # list of bam files
OUT_DIR=$2 # where are our BAM files located?

#   Do the work
parallel indexBAM {} "${OUT_DIR}" :::: "${BAM_LIST}"

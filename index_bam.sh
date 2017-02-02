#!/bin/bash

#PBS -l mem=1gb,nodes=1:ppn=8,walltime=6:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -u
set -o pipefail

#   Define usage message
#   Note: The following paths to arguments will need to be hardcoded \n\
    #   1. [bamFile_list] is a list of bam files to be indexed
    #   2. [bamFile_dir] is where our BAM files are located

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
    samtools index -b "${out_dir}/${sampleName}.bam"
}

#   Export function
export -f indexBAM

#   Arguments provided by user
#   list of bam files
BAM_LIST=/panfs/roc/scratch/liux1299/WBDC_Inversions_Project/seq_handling_parts_ref/WBDC_LianaPop/SAM_Processing/SAMtools/WBDC_LianaPop_Finished_BAM_list.txt
#   where are our BAM files located?
OUT_DIR=/panfs/roc/scratch/liux1299/WBDC_Inversions_Project/seq_handling_parts_ref/WBDC_LianaPop/SAM_Processing/SAMtools/Finished

#   Do the work
parallel indexBAM {} "${OUT_DIR}" :::: "${BAM_LIST}"

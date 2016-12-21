#!/bin/bash

#PBS -l mem=3gb,nodes=1:ppn=16,walltime=24:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -u
set -o pipefail

#   Script written by Chaochih Liu
#   December 21, 2016

#   Define usage message
function usage() {
    echo -e "\
$0: \n\
\n\
This script uses fastq-dump. It takes .sra files and outputs gzipped fastq files. \n\
" >&2
exit 1
}

if [[ $# -lt 3 ]]; then usage; fi

#   Dependencies
module load sratoolkit
module load parallel

#   Where are the .sra files located?
SRA_FILES='/panfs/roc/scratch/pmorrell/WGS_Barley'
#   Where do we want our output to go?
#   As of right now, this part is hardcoded, will eventually fix
OUT_DIR='/home/morrellp/liux1299/scratch/WBDC_Inversions_Project/WGS_raw_fastq'

#   Dump fastq files into out directory
#   The -F option ensures it contains only original sequence name
parallel --dry-run 'fastq-dump --split-files -F --gzip {} --outdir ${OUT_DIR}' ::: ${SRA_FILES}/B1K_04_12.sra
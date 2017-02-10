#!/bin/bash

#PBS -l mem=1gb,nodes=1:ppn=4,walltime=2:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -u
set -o pipefail

#   Define usage message

#   Dependencies
module load bedtools/2.17.0 # Use this version because higher version produces weird output files
module load datamash/1.1.0
module load parallel

#   Function to calculate coverage
#   Concept of command lines borrowed from Li Lei (https://github.com/lilei1/Utilites/blob/master/coverage_bam_cmd.txt)
function calcCoverage() {
    local bam_file="$1"
    local bam_dir="$2"
    local region_file="$3"
    local out_dir="$4"
    #   Sample name
    sampleName=`basename "${bam_file}" .bam`
    #   Generate coverage hist
    bedtools coverage -hist -abam "${bam_dir}/${sampleName}.bam" -b "${region_file}" > ${out_dir}/${sampleName}.coverage.hist.txt
    #   Count number of lines WITHOUT coverage info
    grep "all" "${out_dir}/${sampleName}.coverage.hist.txt" | wc -l > tmp_${sampleName}_noCov.txt
    #   Calculate quantiles and mean for hist files
    head -n -$(cat "${out_dir}/tmp_${sampleName}_noCov.txt") "${out_dir}/${sampleName}.coverage.hist.txt" | \
        datamash --no-strict min 4 q1 4 median 4 mean 4 q3 4 max 4 > ${out_dir}/${sampleName}_coverage_summary.txt
    echo ${sampleName}
    echo "Number of lines without coverage info:" $(cat tmp_${sampleName}_noCov.txt)
    echo "Min   Q1    Median    Mean    Q3    Max"
    cat ${out_dir}/${sampleName}_coverage_summary.txt
}

#   Export function
export -f calcCoverage

#   Arguments provided by user
BAM_FILE_LIST=$1 # list of bam files
BAM_DIR=$2 # Directory our bamfiles are located
REGION_FILE=$3 # Exome capture regions we are using to calculate coverage
OUT_DIR=$4 # Where do our output files go?

#   Array of sample names
#   This code line borrowed from sequence handling (https://github.com/MorrellLAB/sequence_handling/blob/master/Handlers/Coverage_Mapping.sh)
SAMPLE_NAMES=($(xargs --arg-file="${BAM_FILE_LIST}" -I @ --max-args=1 basename @ .bam))
#   Calculate coverage
parallel --jobs 2 --xapply calcCoverage {1} "${BAM_DIR}" "${REGION_FILE}" "${OUT_DIR}" :::: "${BAM_FILE_LIST}" ::: "${SAMPLE_NAMES[@]}"
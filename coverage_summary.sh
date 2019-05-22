#!/bin/bash

set -e
set -u
set -o pipefail

#   Define usage message
function usage() {
    echo -e "\
$0: \n\
\n\
Usage: ./coverage_summary.sh [bam_list] [bam_dir] [region_file] [out_dir] \n\
\n\
NOTE: arguments must be provided in this order. \n\
\n\
where: \n\
1. [bam_list] is a list of bam files to calculate coverage summary \n\
2. [bam_dir] is where our current BAM files are located \n\
3. [region_file] includes regions we are interested in (i.e. BED file) \n\
4. [out_dir] is where we want our files to go \n\
\n\
Required dependencies:
1. Bedtools v2.28.0 or greater
2. datamash v1.1.0 or greater
3. parallel
    " >&2
exit 1
}

if [[ $# -lt 3 ]]; then usage; fi

#   User provided input arguments
#   List of bam files
BAM_LIST=$1
#   Exome capture regions we are using to calculate coverage
REGION_FILE=$2
#   Where do our output files go?
OUT_DIR=$3
#   What is the name of our project?
PROJECT=$4

# Define array to store list of BAM files
declare -a BAM_ARRAY=($(cat "${BAM_LIST}"))

# Check if output directories exist, if not make them
mkdir -p ${OUT_DIR} ${OUT_DIR}/Coverage_Summary ${OUT_DIR}/Coverage_Summary/Histograms

#   Function to calculate coverage
#   Command lines adapted from sequence_handling
#   (https://github.com/MorrellLAB/sequence_handling/blob/master/Handlers/Coverage_Mapping.sh)
function Calc_Coverage() {
    local bam_file="$1"
    local region_file="$2"
    local out_dir="$3"
    local project="$4"
    #   Sample name
    sampleName=$(basename "${bam_file}" .bam)
    #   Generate coverage hist
    bedtools coverage -hist -abam "${bam_file}.bam" -b "${region_file}" > ${out_dir}/Histograms/${sampleName}.hist
    #   Pull out histogram summarizing coverage among "all" features in A
    grep ^all ${out_dir}/Histograms/${sampleName}.hist > ${out_dir}/Histograms/${sampleName}.hist.all.txt

    #   Begin calculating statistics per bp
    #   The minimum is the coverage on the first line of the "all" fields since they're already sorted
    min=$(head -n 1 ${out_dir}/Histograms/${sampleName}.hist.all.txt | awk -F "\t" '{ print $2 }')
    #   The maximum is the coverage on the last line of the "all" fields
    max=$(tail -n 1 ${out_dir}/Histograms/${sampleName}.hist.all.txt | awk -F "\t" '{ print $2 }')
    #   The mean is the sum of (each coverage * the percent of the genome at that coverage)
    mean=$(awk '{ sum += $2*$5 } END { print sum }' ${out_dir}/Histograms/${sampleName}.hist.all.txt)
    #   The mode is the coverage that has the highest percent of the genome at that coverge (excludes zero coverage)
    mode=$(tail -n +2 ${out_dir}/Histograms/${sampleName}.hist.all.txt | sort -grk5,5 | head -1 | awk -F "\t" '{ print $2 }')

    #   The quantiles are a bit tricky...
    #   row_count will count how many rows down the "all" fields we are
    row_count="0"
    #   freq_sum will be the sum of the frequency fields (column 5) from row 0 to row_count
    freq_sum="0"

    #   While freq_sum < 0.25
    while [ $(echo "if (${freq_sum} < 0.25) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        #   freq is the value of the frequency field (column 5) on the row corresponding to row_count
        freq=$(head -n ${row_count} ${out_dir}/Histograms/${sampleName}.hist.all.txt | tail -1 | awk -F "\t" '{print $5}')
        #   Add freq to freq_sum until the while loop exits
        freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done

    #   The first quantile is the coverage on the row at which the cumulative frequency hits 0.25 or greater
    Q1=$(head -n ${row_count} ${out_dir}/Histograms/${sampleName}.hist.all.txt | tail -1 | awk -F "\t" '{print $2}')

    #   Repeat for Q2 (median)
    while [ $(echo "if (${freq_sum} < 0.5) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        freq=$(head -n ${row_count} ${out_dir}/Histograms/${sampleName}.hist.all.txt | tail -1 | awk -F "\t" '{print $5}')
        freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done

    Q2=$(head -n ${row_count} ${out_dir}/Histograms/${sampleName}.hist.all.txt | tail -1 | awk -F "\t" '{print $2}')

    #   Repeat for Q3
    while [ $(echo "if (${freq_sum} < 0.75) 1 else 0" | bc) -eq 1 ]
    do
        ((row_count += 1))
        freq=$(head -n ${row_count} ${out_dir}/Histograms/${sampleName}.hist.all.txt | tail -1 | awk -F "\t" '{print $5}')
        freq_sum=$(echo "${freq_sum} + ${freq}" | bc -l)
    done

    Q3=$(grep "all" "${out_dir}/Histograms/${sampleName}.hist" | head -n ${row_count} | tail -1 | awk -F "\t" '{print $2}')

    #   Append the statistics to the summary file
    echo -e "${sampleName}"'\t'"${min}"'\t'"${Q1}"'\t'"${mode}"'\t'"${Q2}"'\t'"${mean}"'\t'"${Q3}"'\t'"${max}" >> "${out_dir}/${project}_coverage_summary_unfinished.txt"
    #   Put a call to plotCoverage here
}

#   Export function
export -f Calc_Coverage

#   Calculate coverage
parallel --jobs 2 --xapply Calc_Coverage {1} "${BAM_DIR}" "${REGION_FILE}" "${OUT_DIR}" :::: "${BAM_LIST}" ::: "${SAMPLE_NAMES[@]}"

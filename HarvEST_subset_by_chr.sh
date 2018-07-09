#!/bin/bash

set -e
set -o pipefail

#   Define usage mesage
function usage() {
    echo -e "\
$0: \n\
\n\
Usage: ./HarvEST_subset_by_chr.sh [snp_bac] [start_range] [end_range] [out_prefix] [work_dir] \n\
NOTE: arguments must be provided in this order! \n\
\n\
where: \n\
1. [snp_bac] is a file that includes the Query_SNP name and physical position \n\
2. [start_range] is starting range for the number of chromosomes we want to pull down \n\
3. [end_range] is ending chromosome number \n\
(i.e. if chromosomes are labeled chr1H, chr2H, and chr3H and we want all 3 chromosomes, our start range would be 1 and our end range would be 3)
4. [out_prefix] is what prefix will our output filename look like? \n\
5. [work_dir] is where we want to output our files? \n\
" >&2
exit 1
}

if [[ $# -lt 3 ]]; then usage; fi

#   Subset by chr number
function chrSubset() {
    local snp_bac=$1
    local start_range=$2
    local end_range=$3
    local work_dir=$4
    for i in $(seq $start_range $end_range); do grep "chr${i}H" ${snp_bac} > ${work_dir}/tmp_chr${i}H.txt; done
}

export -f chrSubset

#   Add header to subsetted data
function addHeader() {
    local snp_bac=$1
    local start_range=$2
    local end_range=$3
    local out_prefix=$4
    local work_dir=$5
    #   Create header line
    head -n 1 ${snp_bac} > ${work_dir}/header.txt
    echo "Header of SNP_BAC file looks like: \n"
    #   What does our header line look like?
    head ${work_dir}/header.txt
    #   Add header line to beginning of all subsetted data
    for i in $(seq $start_range $end_range); do cat ${work_dir}/header.txt ${work_dir}/tmp_chr${i}H.txt > ${work_dir}/${out_prefix}_Chr${i}H.txt; done
}

export -f addHeader

#   Cleanup temporary files
function cleanUp() {
    local snp_bac=$1
    local start_range=$2
    local end_range=$3
    local work_dir=$4
    #   Remove header file generated
    rm ${work_dir}/header.txt
    #   Remove temporary subsetted data
    for i in $(seq $start_range $end_range); do rm ${work_dir}/tmp_chr${i}H.txt; done
}

export -f cleanUp

#   Arguments provided by user
SNP_BAC=$1 # SNP_BAC.txt pulled from HarvEST
START_RANGE=$2 # Starting chr number
END_RANGE=$3 # Ending chr number
OUT_PREFIX=$4 # output file prefix
WORK_DIR=$5 # where is our working directory?

#   Do the work
#   Subset by chromosome
chrSubset ${SNP_BAC} ${START_RANGE} ${END_RANGE} ${WORK_DIR}
#   Add headers to subsetted data
addHeader ${SNP_BAC} ${START_RANGE} ${END_RANGE} ${OUT_PREFIX} ${WORK_DIR}
#   Remove temporary intermediate files
cleanUp ${SNP_BAC} ${START_RANGE} ${END_RANGE} ${WORK_DIR}

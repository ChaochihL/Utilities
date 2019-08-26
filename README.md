# Utilities

This repository is a collection of miscellaneous scripts that perform various small tasks.

## Contains the following scripts:

| Script | Quick Description |
| ------ | ----------------- |
| `bed_file_summary.R` | Provides mean, five number summary and total length of capture region size of a BED file. |
| `coverage_summary.sh` | Generates Min, Q1, Median, Mean, Q3, and Max coverage summaries using `bedtools` given a list of BAM files and a BED file. |
| `get_reverse_complement` | Directory that contains scripts and documentation to get the reverse complement of sequences. |
| `HarvEST_subset_by_chr.sh` | Takes a HarvEST SNP_BAC.txt file that has PhysPos and Chr_2016 (see `SNP_PhysPosition_Matching.R` script for merging by PhysPos) columns and subsets the data by chromosome. The script returns a new file for every subset group. |
| `HarvEST_to_BED.R` | Takes a HarvEST file that has PhysPos and Chr_2016 (see `SNP_PhysPosition_Matching.R` script for merging by PhysPos) info based on the new barley reference genome and creates a BED file with 4th column containing SNP ID. The output file was intented to be passed to [`bedops`](http://bedops.readthedocs.io/en/v2p4p21/index.html). |
| `LD_data_prep.sh` | Prior to running `LDheatmap.R`, run this script to: 1) Generate sample list containing snps of interest. 2) Make sure genotyping data and sample lists are sorted correctly and duplicates are filtered out. |
| `LDheatmap.R` | This script generates LD heatmaps using genotyping data and physical position information. Please run `LD_data_prep.sh` prior to running this script. |
| `makeGenotypes.R` | This script converts columns in genotype data frame to a format accepted by R package, LDheatmap. |
| `random_sampler.R` | Randomly samples names/values from a given list |
| `renaming_files` | This directory contains file renaming scripts. See readme in the directory for more details. |
| `SNP_PhysPosition_Matching.R` | Takes VCF file (containing physical position, SNP ID, and Chr info) and HarvEST file and adds Physical Position and Chr_2016 columns based on matching SNP names. Note: this script discards any SNPs that do not have matching SNP names. |
| `snpBAC_NA_filtering.R`| Removes all rows with empty cells in Query_SNP column and missing info for PhysPos column. |
| `transpose_data.R`| Transposes data frames |

# Utilities

This repository is a collection of generalized scripts that perform various small tasks.

===

### Contains the following scripts:

- `bed_file_summary.R`: provides mean, five number summary and total length of capture region size of a BED file.

- `HarvEST_to_BED.R`: takes a HarvEST file that has PhysPos and Chr_2016 (see `SNP_PhysPosition_Matching.R` script for merging by PhysPos) info based on the new barley reference genome and creates a BED file with 4th column containing SNP ID. The output file was intented to be passed to [`bedops`](http://bedops.readthedocs.io/en/v2p4p21/index.html).

- `random_sampler.R`: randomly samples names/values from a given list

- `SNP_PhysPosition_Matching.R`: tales VCF file (containing physical position, SNP ID, and Chr info) and HarvEST file and adds Physical Position and Chr_2016 columns based on matching SNP names. Note: this script discards any SNPs that do not have matching SNP names.

- `transpose_data.R`: transposes data frames

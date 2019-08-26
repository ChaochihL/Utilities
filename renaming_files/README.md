# Renaming Files

This directory contains scripts used to rename samples along with some example files and directories for you to test out the scripts.

## Scripts in this directory

| Script | Quick Description |
| ------ | ----------------- |
| `rename_files.py` | Renames files downloaded from NCBI SRA or any file type based on a lookup table that is generated from NCBI SRA Run Selector info table or manually generated. The script has two modes: 1) Dry run, 2) Rename files. See below for more detailed description. |
| `deprecated/file_rename_check.py` | Run this script prior to running [`file_rename.py` on Peter Morrell's GitHub](https://github.com/pmorrell/Utilities/blob/master/file_rename.py) to check that original and new sample names are associated with the correct sample. |

### How does the `rename_files.py` script work?

Below is what an example dry-run looks like.

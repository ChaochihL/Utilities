# Renaming Files

This directory contains scripts used to rename samples along with some example files and directories for you to test out the scripts.

## Scripts in this directory

| Script | Quick Description |
| ------ | ----------------- |
| `rename_files.py` | Renames files downloaded from NCBI SRA or any file type based on a lookup table that is generated from NCBI SRA Run Selector info table or manually generated. The script has two modes: 1) Dry run, 2) Rename files. See below for more detailed description. |
| `deprecated/file_rename_check.py` | Run this script prior to running [`file_rename.py` on Peter Morrell's GitHub](https://github.com/pmorrell/Utilities/blob/master/file_rename.py) to check that original and new sample names are associated with the correct sample. |

An example directory with filenames to test the script is located in the `toy_rename_bam` directory. I have included a few example lookup tables (two with SRA IDs and one with BAM accession names) located in the `example_lookup_tables` directory. Note: if the line starts with `#`, the script will skip it when reading in the files, so you can use the `#` to start column header lines.

### How does the `rename_files.py` script work?

Below is an example dry-run for a directory with toy BAM files (note: BAM files are empty files since we only need the filename to do a test rename run). The BAM file names are in the directory called `toy_rename_bam`. The example lookup table we will use is `example_lookup_tables/toy_rename_bam_lookup_table.txt`.

```bash
# Running from within the renaming_files directory
# Dry-run to see what will be renamed
# Note: we are only showing the first few lines here
./rename_files.py example_lookup_tables/toy_rename_bam_lookup_table.txt toy_rename_bam/ --dry-run
Dry-run, print old name and new name. Please run with --rename option to do the actual renaming.
Old_Name	 New_Name
CIho_00497_S24.bai CIho_00497.bai
CIho_00497_S24.bam CIho_00497.bam
CIho_04083_S2.bai CIho_04083.bai
CIho_04083_S2.bam CIho_04083.bam
CIho_01468_S1.bai CIho_01468.bai
CIho_01468_S1.bam CIho_01468.bam
CIho_04088_S4.bai CIho_04088.bai
CIho_04088_S4.bam CIho_04088.bam
CIho_14059_S3.bai CIho_14059.bai
CIho_14059_S3.bam CIho_14059.bam
PI_038320_S5.bai PI_038320.bai
PI_038320_S5.bam PI_038320.bam
PI_061589_S6.bai PI_061589.bai
PI_061589_S6.bam PI_061589.bam
PI_064004_S23.bai PI_064004.bai
PI_064004_S23.bam PI_064004.bam
.
.
.
```

The dry-run mode serves as a check before doing the actual renaming. Here, we use a path relative to our current working directory. We can also specify absolute filepaths to the lookup tables and the directory containing files to rename.

To do the actual renaming, add the `--rename` option.

# Aligned vs FASTQ Data Verification

This directory contains a set of scripts to check aligned (SAM/BAM) files are correct and did not get partially overwritten during the read mapping process. These scripts combined compares the Illumina sequence identifiers from the aligned file to Illumina sequence identifiers in the raw FASTQ files. If there are mismatches in sequence identifiers for the same accessions, the scripts check all other FASTQ files to find where the sequence identifier came from.

**Important note:** The scripts currently assume the input FASTQ files are gzipped, we are working with paired end reads, and that BWA was used for read mapping.

---

### Usage

Fill out filepaths to input arguments in the `driver_check_aligned_vs_fastq.job` script. This script then calls the `bam_vs_fastq_check.sh` script, which calls the `check_seq_id.py` script.

There are two modes to run in: 1) `CHECK_SEQIDS` and 2) `SEQID_ORIGIN`.
'CHECK_SEQIDS' outputs a summary of proportion mismatch for each accession and files that tell you which sequence identifers matched/mismatched between SAM/BAM file and fastq file.
'SEQID_ORIGIN' does the same as 'CHECK_SEQID' but adds an additional search to find the origin of the mismatched sequence identifiers. This can be informative if there are any file name mixups.

After filling out the driver script, run the following:

```bash
# Submit as job array to PBS HPC system
# Where "n" is: the maximum number of samples in ACC_LIST - 1
# Ex: If you have 24 accessions in ACC_LIST your range would be 0-23
qsub -t 0-n driver_check_aligned_vs_fastq.job
```

### Dependencies
- Python 3 v3.7 or later (specific libraries below)
    - argparse
    - sys
    - os
    - re
    - Biopython
    - gzip
- Samtools v1.8.0 or later


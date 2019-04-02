"""This script compares Illumina sequence identifiers in aligned (SAM or BAM)
files to Illumina sequence identifiers in FASTQ files. Script outputs a table
containing a column of sequence identifiers and info on if they come from the
correct sample.

Usage:
python3 check_seq_id.py [accession] [aligned_header] [aligned_seqIDs] [fastq_R1] [fastq_R2]

Where:
x) [accession] is a single accession name.
x) [aligned_header]
x) [aligned_seqIDs] is the full filepath to a file containing Illumina
    sequence identifiers from the SAM/BAM file.
x) [fastq_R1] list of full filepaths for forward gzipped fastq files.
x) [fastq_R2] list of full filepaths for reverse gzipped fastq files.

Note: Script assumes user is running Python v3.7 or later due to the assumption
that dictionaries maintain insertion order. Script also assumes BWA was used
for read mapping. This matters for the function that parses the aligned file
header.
"""

import sys
import os
from Bio import SeqIO
import gzip


def read_fastq(fastq_file):
    """Parse the gzipped FASTQ file and keep only the sequence identifiers
    in a dictionary."""
    fastq_dict = {}
    with gzip.open(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            fastq_dict[record.name] = [record.id, record.description]
    return (fastq_dict)


def read_aligned_header(aligned_header_file):
    """Parse SAM/BAM file header lines and extract @PG header line. @PG header line
    contains the BWA command ran and fastq files that were aligned. This gives
    us info about which sample was actually processed."""
    with open(aligned_header_file, "r") as handle:
        for line in handle:
            if "@PG" and "bwa" in line:
                tmp_pg = line.split()
    # Pull out fastq files bwa actually aligned
    fastq_filenames = []
    for i in range(1, len(tmp_pg)):
        if fastq_suffix in tmp_pg[i]:
            fastq_filenames.append(os.path.basename(tmp_pg[i]))
    return (fastq_filenames)


def read_aligned_seqids(aligned_seqid_file):
    """Parse a file containing Illumina sequence identifiers from a
    SAM/BAM file. This file is created in the bam_vs_fastq_check.sh
    script, which calls on this script."""
    aligned_seqids = []
    with open(aligned_seqid_file, "r") as handle:
        for line in handle:
            aligned_seqids.append(line.strip('\n'))
    return(aligned_seqids)


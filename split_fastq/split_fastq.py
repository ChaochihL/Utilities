#!/usr/bin/env python3

"""Script that splits fastq files into N parts (specified by the user).

Usage:
./split_fastq.py [fastq_file_path] [num_parts] [total_num_seqs] [out_dir_fp]

Where:
1) [fastq_file_path] is the full filepath to the input FASTQ file. FASTQ file
        can be uncompressed or bgzipped. IMPORTANT: If files are NOT compressed
        using bgzip, the script will return an error.

2) [num_parts] How many pieces do we want to split our file into? (integer)

3) [total_num_seqs] How many total sequences do we have? You can get this
        number by running FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        and looking at the value for "Total Sequences" in the HTML report.

4) [out_dir_fp] is the full filepath to where we will store the split fastq
        files.
"""

import sys
import os
import math
import gzip

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(0)


def prep_file_prefix(fastq_filepath):
    """Generate prefix for split output files."""
    # Strip file extension from input fastq file
    if ".gz" in fastq_filepath:
        # Strip .gz part of file extension
        temp_name = os.path.splitext(os.path.basename(fastq_filepath))[0]
        # Strip .fastq or .fq part of file extension
        file_prefix = os.path.splitext(os.path.basename(temp_name))[0]
    else:
        # Assume file extension is either .fastq or .fq
        file_prefix = os.path.splitext(os.path.basename(fastq_filepath))[0]
    return file_prefix


def prep_out_fnames(file_prefix, num_parts, num_seqs_per_part):
    """Generate a list of output filenames and the maximum number of sequences
    in each output file."""
    out_files = []
    curr_num_seqs = 0
    for i in range(num_parts):
        # Want parts to start with part 1 instead of part 0
        curr_part = i + 1
        curr_filename = file_prefix + "_part" + str(curr_part) + ".fastq"
        starting_num_seqs = curr_num_seqs
        curr_num_seqs = curr_num_seqs + num_seqs_per_part
        out_files.append([curr_filename, starting_num_seqs, curr_num_seqs])
    # Subtract one from upper bound to make the counting easier
    for i, elem in enumerate(out_files):
        out_files[i][2] = out_files[i][2] - 1
    return out_files


def split_fastq_file(f, out_files, out_dir):
    """Split our fastq file into N parts and save to file. We also minimize
    the number of times the output part files are opened/closed to reduce
    overhead. Each part file is only opened once and closed once."""
    # Split file and write to parts
    findex = 0  # Start at file index 0
    curr_count = 0
    # Iterate line by line
    for i in f:
        if (curr_count >= out_files[findex][1] and curr_count <= out_files[findex][2]):
            # print("Current file index is:", findex)
            # Output sequence to current part
            curr_out_fp = out_dir + "/" + out_files[findex][0]
            # Open file (one time) for writing if file doesn't exist
            if not os.path.exists(curr_out_fp):
                # Create and open file
                print("Opening file:", curr_out_fp)
                cof = open(curr_out_fp, 'a+')
            # Write line to file
            cof.write(i)
            # Increment by one each sequence (starts with "@")
            if i.startswith("@"):
                curr_count += 1
            # print("after writing sequence count:", curr_count)
        elif curr_count > out_files[findex][2] and not i.startswith("@"):
            # print("Current file index is:", findex)
            # Write line to file
            cof.write(i)
            # print("after writing sequence count:", curr_count)
        elif curr_count > out_files[findex][2] and i.startswith("@"):
            # Check if we are processing our last sequence
            if findex == len(out_files) - 1:
                # Don't increment file index anymore
                print("We are at the last sequence.")
            else:
                # Close previous part file before incrementing
                print("Closing file:", curr_out_fp)
                cof.close()
                # Move on to next output file
                findex += 1
            # print("Current file index is:", findex)
            # Output sequence to current part
            curr_out_fp = out_dir + "/" + out_files[findex][0]
            # print(curr_out_fp)
            # Open file (one time) for writing
            if not os.path.exists(curr_out_fp):
                # Create and open file
                print("Opening file:", curr_out_fp)
                cof = open(curr_out_fp, 'a+')
            # Write line to file
            cof.write(i)
            # Increment by one each sequence (starts with "@")
            if i.startswith("@"):
                curr_count += 1
            # print("after writing sequence count:", curr_count)
    # Finished iterating through all sequences, close last part file
    print("Closing file:", curr_out_fp)
    cof.close()


def main(fastq_file_path, num_parts, total_num_seqs, out_dir_fp):
    """Driver function."""
    # Prepare filepaths
    fastq_fp = os.path.expanduser(fastq_file_path)
    # Strip trailing slash so output path doesn't get messed up
    out_dir = os.path.expanduser(out_dir_fp).rstrip("/")
    # Prepare file prefix for split output files
    file_prefix = prep_file_prefix(fastq_fp)
    # Identify the number of sequences to output in each part
    num_seqs_per_part = math.ceil(int(total_num_seqs) / int(num_parts))
    # Generate list of output filenames
    out_files = prep_out_fnames(file_prefix, int(num_parts), num_seqs_per_part)
    # Prior to splitting file, clean up any existing parts (since we
    # are appending later)
    for l in out_files:
        curr_out_fp = out_dir + "/" + l[0]
        if os.path.exists(curr_out_fp):
            print("More than 1 part file exists, remove to start from clean slate:", curr_out_fp)
            os.remove(curr_out_fp)
        else:
            print(l[0], "file doesn't exist, proceeding to split file...")
    # Split fastq file into parts and save to file
    # Read in file line by line
    if "gz" in fastq_fp:
        with gzip.open(fastq_fp, 'rt') as f:
            split_fastq_file(f, out_files, out_dir)
    else:
        with open(fastq_fp, 'r') as f:
            split_fastq_file(f, out_files, out_dir)
    print("Done splitting fastq files.")


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

#!/usr/bin/env python3
# Chaochih Liu - Falcon Heights, MN June 17, 2019
"""This script renames SRA/fastq files obtained from NCBI SRA.
This script takes in a tab-delimited lookup table with two columns:
    1) Run/Exp/SRA ID
    2) Corresponding accession name

The lookup table can be easily created from NCBI SRA Run Selector and downloading the run info table. Search for NCBI SRA and enter in the BioProject ID (ex: PRJNA399170), then send the results to run selector.

To pull up usage message, run: ./rename_sra_fastq.py --help
"""

import os
import sys
import argparse


def parse_args():
    """Set up argument parser to parse command line options."""
    parser = argparse.ArgumentParser(
        description=('Rename SRA/fastq files obtained from NCBI SRA.'
                     'This script takes in a tab-delimited lookup table with two columns:'
                     '1) Run/Exp/SRA ID, 2) Corresponding accession name. '
                     'The lookup table can be easily created from NCBI SRA Run Selector and downloading '
                     'the run info table. Search for NCBI SRA and enter in the BioProject ID '
                     '(ex: PRJNA399170), then send the results to run selector. Note: Script currently '
                     'works on entire directories and every file prefix needs to be in the lookup table.'),
        add_help=True
    )
    # Define required arguments
    parser.add_argument(
        'lookup_table',
        metavar='lookup_table',
        help='Full filepath to lookup table with two columns: 1) Run/Exp/SRA ID, 2) Accession name.'
    )
    parser.add_argument(
        'file_directory',
        metavar='file_dir',
        help='Full filepath to directory containing files to be renamed.'
    )
    # Define optional arguments
    # So that user can either run a dry-run mode or actual renaming mode
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        '--dry-run',
        action='store_true',
        help='Check before actual run. Prints old name and new name.'
    )
    mode.add_argument(
        '--rename',
        action='store_true',
        help='Rename the samples based on lookup table naming pairs.'
    )
    a = parser.parse_args()
    return a


def set_mode(args):
    """Set the mode for which check will be run. Options are:
    1) --dry-run: A quick check before actual run. Prints old name and new name.
    2) --rename: does the actual file renaming."""
    arg_dict = vars(args)
    if arg_dict['dry_run']:
        print('Dry-run, print old name and new name. Please run with --rename option to do the actual renaming.')
        return ('DRY_RUN')
    elif arg_dict['rename']:
        print('Renaming files...')
        return ('RENAME')
    else:
        print('Please specify either --dry-run or --rename options.')
        return(None)


def read_table(lookup_table):
    """Parse the lookup table and store in a dictionary."""
    lt_dict = {}
    with open(lookup_table, "rt") as handle:
        for line in handle:
            # Skip the header lines, if there are any
            if line.startswith('#'):
                continue
            else:
                # Strip new line and split by tab delimiter
                tmp = line.strip().split('\t')
                lt_dict[tmp[0]] = tmp[1]
    return lt_dict


def rename_file(lt_dict, files):
    """Create paired rename list."""
    rename_dict = {}
    for k in lt_dict:
        for i in files:
            if k in i:
                # Create dictionary with old name and new name
                rename_dict[i] = i.replace(k, lt_dict[k])
    return rename_dict


def main():
    """Driver function."""
    args = parse_args()
    f = set_mode(args)

    # Read in lookup table
    lookup_dict = read_table(args.lookup_table)
    # Create list of files
    files = sorted(os.listdir(args.file_directory))
    # Prepare to rename file
    rename_dict = rename_file(lookup_dict, files)

    # Parser user provided arguments and run requested options
    if f == 'DRY_RUN':
        # Dry-run, print old name and new name only
        print("Old_Name\t", "New_Name")
        for k in rename_dict:
            print(k, rename_dict[k])
    elif f == 'RENAME':
        os.chdir(args.file_directory)
        for k in rename_dict:
            os.rename(k, rename_dict[k])
    else:
        print('Please specify --rename option to do actual renaming.')
        exit(1)
    return


main() # Run the program

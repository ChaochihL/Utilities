#!/usr/bin/env python3
"""Compares Illumina sequence identifiers in aligned (SAM or BAM) files to
Illumina sequence identifiers in FASTQ files. Script outputs a table containing
a column of sequence identifiers and info on if they come from the correct
sample.

For usage message, run: ./check_seq_id.py --help

Note: Script assumes BWA was used for read mapping. This matters for the
function that parses the aligned file header.
"""

import sys
import os
import re
import argparse
from Bio import SeqIO
import gzip


def parse_args():
    """Set up argument parser to parse command line options."""
    parser = argparse.ArgumentParser(
        description="Check SAM/BAM vs FASTQ file sequence identifiers.",
        add_help=True
    )
    # Define required arguments
    parser.add_argument(
        'accession',
        metavar='acc',
        help='A single accession name from aligned file.'
    )
    parser.add_argument(
        'aligned_header_file',
        metavar='aligned_header',
        help=('Full filepath to a file containing SAM/BAM header lines '
              'for accession we are processing.')
    )
    parser.add_argument(
        'aligned_seqIDs_file',
        metavar='aligned_seqIDs',
        help=('File containing Illumina sequence sequence identifiers '
              'from the SAM/BAM file.')
    )
    parser.add_argument(
        'fastq_R1_file',
        metavar='fastq_R1',
        help='Forward gzipped FASTQ file corresponding to accession'
    )
    parser.add_argument(
        'fastq_R2_file',
        metavar='fastq_R2',
        help='Reverse gzipped FASTQ file corresponding to accession.'
    )
    parser.add_argument(
        'fastq_list_fp',
        metavar='fastq_list_fp',
        help='List of full filepaths for all FASTQ files.'
    )
    parser.add_argument(
        'fastq_suffix',
        metavar='fastq_suffix',
        help='Suffix that matches FASTQ file suffix. Ex: .fastq.gz'
    )
    parser.add_argument(
        'out_dir',
        metavar='out_dir',
        help='Full filepath to output directory where files will be saved.'
    )
    # Define optional arguments
    # So that user can only select one or the other mode of check
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        '--check-seqids',
        action='store_true',
        help=('Checks for mismatches between SAM/BAM and FASTQ sequence IDs. '
              'Only tells you if there is a mismatch or not.')
    )
    mode.add_argument(
        '--seqid-origin',
        action='store_true',
        help=('Searches through list of fastq files to find origin of '
              'mismatched sequence IDs. Output tells you if there is a '
              'mismatch or not and where the mismatched sequence IDs '
              'originated. Important: only include this argument if '
              'you want this extra search to run')
    )
    a = parser.parse_args()
    return (a)


def set_mode(args):
    """Set the mode for which check will be run. Options are:
    1) check-seqids: only does a check and tell you if there are mismatches.
    2) find-seqid-origin: does a check and finds the origin of mismatched
    seqids."""
    arg_dict = vars(args)
    # print(arg_dict)
    if arg_dict['check_seqids']:
        print('Running --check-seqids.')
        return ('CHECK_SEQIDS')
    elif arg_dict['seqid_origin']:
        print('Running --seqid-origin.')
        return ('SEQID_ORIGIN')
    else:
        print('No options specified, please specify')
        return(None)


def read_fastq(fastq_file):
    """Parse the gzipped FASTQ file and keep only the sequence identifiers
    in a dictionary."""
    fastq_dict = {}
    with gzip.open(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            fastq_dict[record.name] = [record.id, record.description]
    return (fastq_dict)


def read_aligned_header(aligned_header_file, fastq_suffix):
    """Parse SAM/BAM file header lines and extract @PG header line. @PG header line
    contains the BWA command ran and fastq files that were aligned. This gives
    us info about which sample was actually processed."""
    with open(aligned_header_file, "r") as handle:
        for line in handle:
            if "@PG" and "bwa" in line:
                tmp_pg = line.split()
    # Pull out fastq files bwa actually aligned
    fastq_filenames = []
    for i in range(0, len(tmp_pg)):
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
    # For paired end reads, there will be duplicate seqids
    # Get unique seqids
    aligned_seqids_uniq = list(set(aligned_seqids))
    return (aligned_seqids_uniq)


def check_acc_names(acc, aligned_header_fastq):
    """Check if the accession name matches fastq file that was aligned
    (shown in SAM/BAM header lines)."""
    if acc in aligned_header_fastq[0]:
        acc_info = [acc, aligned_header_fastq[0], aligned_header_fastq[1],
                    "match"]
    else:
        acc_info = [acc, aligned_header_fastq[0], aligned_header_fastq[1],
                    "mismatch"]
    return (acc_info)


def check_seqids(acc, a_seqids, fastq_r1_dict, fastq_r2_dict):
    """For accession we are currently processing, compare aligned seqIDs
    to fastq file for same accession. Output will have 3 columns:
    1) seqID, 2) match/mismatch, and 3) accession"""
    # Check aligned seqids vs raw fastq r1 seqid
    seqids_r1 = {}
    # Keep track of number of mismatches found to get proportion mismatches
    count_r1 = 0
    for s in range(0, len(a_seqids)):
        if a_seqids[s] in fastq_r1_dict.keys():
            seqids_r1[a_seqids[s]] = [a_seqids[s], "match", acc]
        else:
            seqids_r1[a_seqids[s]] = [a_seqids[s], "mismatch", "un"]
            count_r1 += 1
    # Check aligned seqids vs raw fastq r2 seqid
    seqids_r2 = {}
    # Keep track of number of mismatches found to get proportion mismatches
    count_r2 = 0
    for s in range(0, len(a_seqids)):
        if a_seqids[s] in fastq_r2_dict.keys():
            seqids_r2[a_seqids[s]] = [a_seqids[s], "match", acc]
        else:
            seqids_r2[a_seqids[s]] = [a_seqids[s], "mismatch", "un"]
            count_r2 += 1
    # Calculate proportion of mismatches for accession
    prop_mismatch_r1 = count_r1/len(seqids_r1)
    prop_mismatch_r2 = count_r2/len(seqids_r2)
    prop_mismatch = [acc, str(prop_mismatch_r1), str(prop_mismatch_r2)]
    return (seqids_r1, seqids_r2, prop_mismatch)


def prep_fastq_lists(fastq_list_fp):
    """Assumes forward reads are denoted by 'R1' and reverse reads
    are denoted by 'R2'. Please modify this function if using other
    naming scheme for forward and reverse reads."""
    fastq_list = []
    with open(fastq_list_fp, "r") as f:
        for line in f:
            fastq_list.append(line.strip('\n'))
    # Create forward and reverse fastq list
    # R1
    fastq_r1_list = []
    for f in fastq_list:
        if "R1" in f:
            fastq_r1_list.append(f)
    # R2
    fastq_r2_list = []
    for f in fastq_list:
        if "R2" in f:
            fastq_r2_list.append(f)
    return (fastq_r1_list, fastq_r2_list)


def find_unseqid_origin(acc_info, seqids, fastq_list):
    """For each aligned seqID, search through all fastq files in list
    and find where the aligned seqID came from. Once found, change the "un"
    value to the fastq filename where aligned seqID is found in fastq file"""
    # To reduce total number of searches performed, we will assume that
    # the fastq file aligned (shown in aligned file header) is likely
    # where the mismatched seqids came from. We will start our search there
    # by first pulling out that fastq accession name.
    aligned_fq = acc_info[1].split('_')
    fq_intersect = []
    for f in fastq_list:
        tmp_fq = os.path.basename(f).split('_')
        tmp_set = frozenset(tmp_fq)
        # Checking for membership in a set is roughly O(1) compared to O(n)
        # for list.
        intersection = [x for x in aligned_fq if x in tmp_set]
        # Note: this condition assumes that your accession names are separated
        # by an underscore (i.e., WBDC_336). If your accession name is
        # something like WBDC336, you may want to change the '> 1' to '> 0'.
        if len(intersection) > 1:
            fq_intersect.append('_'.join(intersection))
        else:
            continue
    print("FASTQ file intersection: ", ','.join(fq_intersect))

    for key in seqids:
        # Only search for aligned seqIDs that have unknown origin
        if seqids[key][2] == "un":
            # Use to stop searching through remaining FASTQ files after
            # we have found where the aligned sequence identifier came from
            breaker = "not_found"

            # Search through mismatched sample from aligned header first
            start_fq = []
            # Pull out fastq filepath to start search with
            for fp in fastq_list:
                if len(fq_intersect) == 0:
                    # This part assumes underscore separates accession names
                    # i.e., WBDC_028
                    tmp = os.path.basename(fp)
                    if re.search(aligned_fq[1], tmp):
                        start_fq.append(fp)
                elif len(fq_intersect) > 0 and fq_intersect[0] in fp:
                    start_fq.append(fp)
                else:
                    print("No starting fastq file.")
            with gzip.open(start_fq[0], "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    if seqids[key][2] == "un" and key == record.id:
                        # Extract filename for where aligned seqID
                        # matches seqID in fastq file
                        filename = os.path.basename(start_fq[0])
                        # Replace "un" with filename
                        seqids[key][2] = filename
                        # Change to "found_match" so we stop searching through
                        # remaining fastq files
                        breaker = "found_match"
                        print("Found match, sample is: ", filename)
                        break
            # Check remaining fastq files in list
            for f in range(0, len(fastq_list)):
                # Skip processing start_fq sample and self sample
                if (
                    breaker == "not_found" and
                    start_fq[0] not in fastq_list[f] and
                    acc_info[0] not in fastq_list[f]
                ):
                    print(acc_info[0], key, "not found, continue search")
                    with gzip.open(fastq_list[f], "rt") as handle:
                        for record in SeqIO.parse(handle, "fastq"):
                            if seqids[key][2] == "un" and key == record.id:
                                # Extract filename for where aligned seqID
                                # matches seqID in fastq file
                                filename = os.path.basename(fastq_list[f])
                                # Replace "un" with filename
                                seqids[key][2] = filename
                                # Change to "found_match" so we stop searching
                                # through remaining fastq files
                                breaker = "found_match"
                                print("Found match, sample is: ", filename)
                                break
                else:
                    continue
        else:
            continue
    return (seqids)


def save_to_file(accession_info, seqids, outfile_name):
    """Save data to output file."""
    header_line = ["#SeqID", "Aligned_vs_Fastq_SeqID", "Fastq_Origin"]
    with open(outfile_name, 'a') as out:
        # Define info lines at top of file
        out.write("#Info: " + ', '.join(accession_info) + '\n')
        out.write("#Info line describes in the following order: 1) Aligned "
                  "accession, 2) and 3) FASTQ files aligned, "
                  "4) match/mismatch between accession and FASTQ file "
                  "aligned." + '\n')
        out.write('\t'.join(header_line) + '\n')
        # Save seqids check to output file
        for key in seqids:
            out.write('\t'.join(seqids[key]) + '\n')
    return


def save_prop_mismatch(prop_mismatch, outfile_name):
    """Save proportion mismatch for accession to a summary table.
    This function processes one accession at a time, however specifying
    same filename will append new accession summaries to existing file."""
    header_line = ["#Accession", "R1_prop_mismatch", "R2_prop_mismatch"]
    # Check if output file exists
    tmp_exists = os.path.isfile(outfile_name)
    if tmp_exists:
        # Check if header line exists in file
        tmp_file = open(outfile_name, 'r')
        # Check first line
        if "#Accession" not in tmp_file.readline():
            with open(outfile_name, 'a') as out:
                # Add header line
                out.write('\t'.join(header_line) + '\n')
                # Then add proportion mismatch
                out.write('\t'.join(prop_mismatch) + '\n')
        else:
            with open(outfile_name, 'a') as out:
                out.write('\t'.join(prop_mismatch) + '\n')
    else:
        with open(outfile_name, 'a') as out:
            # Add header line
            out.write('\t'.join(header_line) + '\n')
            # Then add proportion mismatch
            out.write('\t'.join(prop_mismatch) + '\n')
    return


def driver_check_seqids(accession, aligned_header_file, aligned_seqIDs_file,
                        fastq_R1_file, fastq_R2_file, fastq_suffix):
    """Driver function that reads in all files necessary and checks for
    mismatches between SAM/BAM file sequence identifiers and FASTQ file
    sequence identifiers. Outputs from this function will only tell you
    if there is a mismatch and which sequence identifiers have a mismatch.
    It will not tell you where the mismatched sequence identifiers
    originated."""
    # Expand user filepaths
    aligned_header_fp = os.path.expanduser(aligned_header_file)
    aligned_seqIDs_fp = os.path.expanduser(aligned_seqIDs_file)
    fastq_r1_fp = os.path.expanduser(fastq_R1_file)
    fastq_r2_fp = os.path.expanduser(fastq_R2_file)

    # Read in files
    fastq_r1_dict = read_fastq(fastq_r1_fp)
    fastq_r2_dict = read_fastq(fastq_r2_fp)
    aligned_header = read_aligned_header(aligned_header_fp, fastq_suffix)
    aligned_seqids = read_aligned_seqids(aligned_seqIDs_fp)

    # Check accession vs fastq file that was actually processed
    acc_info = check_acc_names(accession, aligned_header)
    # Compare aligned sequence identifiers to those in raw fastq file
    # for accession we are currently working with
    seqids_r1, seqids_r2, prop_mismatch = check_seqids(
        accession,
        aligned_seqids,
        fastq_r1_dict,
        fastq_r2_dict
    )
    return(acc_info, seqids_r1, seqids_r2, prop_mismatch)


def driver_find_seqid_origin(fastq_list_fp, seqids_r1, seqids_r2, acc_info):
    """Driver function that takes the output from driver_check_seqids()
    and performs a search through a list of fastq files to identify
    origin of mismatched sequence identifiers."""
    # Prep fastq lists
    fastq_r1_list, fastq_r2_list = prep_fastq_lists(fastq_list_fp)
    # For aligned seqIDs where it is not in the fastq file for the
    # sample we are currently processing, search through all other
    # fastqs to find where aligned seqID came from
    print("Processing R1, find seqid origin")
    seqids_r1_out = find_unseqid_origin(acc_info, seqids_r1, fastq_r1_list)
    print("Processing R2, find seqid origin")
    seqids_r2_out = find_unseqid_origin(acc_info, seqids_r2, fastq_r2_list)
    return(seqids_r1_out, seqids_r2_out)


def main():
    """Main function."""
    args = parse_args()
    print("Processing accession:", args.accession)
    f = set_mode(args)

    # Parse user provided arguments and run requested checks
    if f == 'CHECK_SEQIDS':
        # Check seqid
        acc_info, seqids_r1, seqids_r2, prop_mismatch = driver_check_seqids(
            args.accession,
            args.aligned_header_file,
            args.aligned_seqIDs_file,
            args.fastq_R1_file,
            args.fastq_R2_file,
            args.fastq_suffix
        )
        # Save proportion mismatches summary to output file
        summary_outname = (args.out_dir + '/' +
                           "temp_" + args.accession +
                           "_prop_mismatch.txt")
        save_prop_mismatch(prop_mismatch, summary_outname)
        # Save seqid checks to output file
        r1_outname = (args.out_dir + '/' + args.accession +
                      "_R1_seq_id_check.txt")
        r2_outname = (args.out_dir + '/' + args.accession +
                      "_R2_seq_id_check.txt")
        # R1
        save_to_file(acc_info, seqids_r1, r1_outname)
        # R2
        save_to_file(acc_info, seqids_r2, r2_outname)
    elif f == 'SEQID_ORIGIN':
        # Check seqid
        acc_info, seqids_r1, seqids_r2, prop_mismatch = driver_check_seqids(
            args.accession,
            args.aligned_header_file,
            args.aligned_seqIDs_file,
            args.fastq_R1_file,
            args.fastq_R2_file,
            args.fastq_suffix
        )
        # Save proportion mismatches summary to output file
        summary_outname = (args.out_dir + '/' +
                           "temp_" + args.accession +
                           "_prop_mismatch.txt")
        save_prop_mismatch(prop_mismatch, summary_outname)
        # Find seqid origin
        seqids_r1_out, seqids_r2_out = driver_find_seqid_origin(
            args.fastq_list_fp,
            seqids_r1,
            seqids_r2,
            acc_info
        )
        # Save seqid checks to output file
        r1_outname = (args.out_dir + '/' + args.accession +
                      "_R1_seq_id_check_origin.txt")
        r2_outname = (args.out_dir + '/' + args.accession +
                      "_R2_seq_id_check_origin.txt")
        # R1
        save_to_file(acc_info, seqids_r1_out, r1_outname)
        # R2
        save_to_file(acc_info, seqids_r2_out, r2_outname)
    else:
        print('Please specify which check to perform:'
              '--check-seqids or --seqid-origin')
        exit(1)
    return

main()

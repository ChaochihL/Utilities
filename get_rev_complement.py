#!/usr/bin/env python3

#   Chaochih Liu - Oct 31, 2017

"""This script takes in an adapters.fa file and adds the reverse complement sequence for each adapter.
Requires biopython.

Usage: ./get_rev_complement [adapters_file]

Where:
    [adapters_file] is a file that has adapter names starting with ">" (i.e. >H516_index2_i5)
    and sequences starting on new line following

    Example format:

    >H516_index2_i5
    ACTCTAGG
"""

import sys
from Bio.Seq import Seq

def read_adapters_file(adapters_file):
    """Function that reads in adapters.fa file."""
    #   Use list comprehension to read in file line by line
    #   Removing trailing new line
    a_dat = [line.strip('\n') for line in open(adapters_file, 'r').readlines()]
    return a_dat


def main(file1):
    """Function that runs the program by reading in file, creating dictionary from adapter name and sequences, and getting reverse complement sequences."""
    #   Read in adapters.fa file from command line
    adapters = read_adapters_file(file1)

    #   Create dictionary with adapter name as key and adapter sequence as value
    adapter_dict = {}
    for i in range(0, len(adapters)):
        if adapters[i].startswith(">"):
            #   If adapter starts with ">", this is the adapter name (key)
            #   The line immediately proceeding it, is the value (hence the i+1)
            adapter_dict[adapters[i]] = adapters[i+1]


    for key, value in adapter_dict.items():
        #   Create a sequence object
        tmp_seq = Seq(adapter_dict[key])
        #   Get reverse complement
        rev_comp = tmp_seq.reverse_complement()
        #   Print original adapter name and adapter sequence
        #   followed by adapter reverse complement name and reverse complement sequence
        print(key, value, key + "_reverse_complement", rev_comp, sep = '\n')


main(sys.argv[1]) # Run the program

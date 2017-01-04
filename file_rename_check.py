#!/usr/bin/env python3

#   Script written by Chaochih Liu
#   November 29, 2016

"""
A script to check that original and new sample names are the same samples.

Meant as a check prior to running Peter Morrell's file_rename.py script.
"""

import re
import sys
import argparse

Usage="""
This script reads in a file with the following columns:
1. Accession name
2. SRA/ENA Accession code
3. Raw data sample names

The output file will consist of two columns:
1. Original sample names
2. New sample names for original data

Usage:
    file_rename_check.py [sample_info.txt]
"""

#   Print usage message if no arguments are provided
if not sys.argv[1:]:
    print(Usage)
    exit(1)


#   Input file should contain columns specified under Usage message
names = sys.argv[1]
names = 'toy_IPK_spont.txt' # test


#   Read in file
my_query = [] # This is our search key
accession = [] # This is the original data
#with open(sys.argv[1], 'r') as accn:
with open(names, 'r') as accn: # test, use line above when done testing
    #   Skip header line
    for line in accn.readlines()[1:]:
        tmp = line.split("\t")
        accession.append((tmp[0], tmp[1], tmp[2]))
        my_query.append((tmp[1]))


#   Check for matching Accession_Code and Original_Name
no_match = []
for i in my_query:
    key = re.compile(i) # create search string
    for a in accession:
        key.findall(a)


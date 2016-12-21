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

if not sys.argv[1:]:
    print(Usage)
    exit(1)



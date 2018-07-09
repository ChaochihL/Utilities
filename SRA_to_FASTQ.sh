#!/bin/bash

set -e
set -u
set -o pipefail

#   Script written by Chaochih Liu
#   December 21, 2016

USAGE="Usage:
$0 [Option] [SRA file list] [ -d DIR ] [ -h ]

Description:
    This script uses the SRA Toolkit to convert .sra files to gzipped fastq files in parallel.
    It takes in a list of full filepaths to the .sra files. If DIR is specified, all output
    files will be put there (script creates DIR if it doesn't exist). If DIR is not
    supplied, then the current directory is used.

Dependencies:
    1) NCBI SRA Toolkit
    2) GNU parallel

Pass -h to see available options.
"

HELP="Usage:
$0 [Option] [SRA file list] [ -d DIR ] [ -h ]

Available options:
Required:
    -l      Provided [SRA file list]
Optional:
    -d DIR  Output all .SRA files into DIR. Provide full filepath to DIR.
Switches:
    -h      Show this message and exit.
"

#   If there are no arguments passed to the script, drop the usage message and exit
if [ $# == 0 ]
    then echo "$USAGE"
    exit 1
fi

#   Parse the options
DIR=$(pwd)
while [[ $# > 0 ]]
    do
        FLAG="$1"
        case $FLAG in
            -l)
            SRA_LIST="$2"
            shift
            ;;
            -d)
            DIR="$2"
            shift
            ;;
            -h)
            echo "$HELP"
            exit 2
            ;;
            *)
            echo "$USAGE"
            exit 1
            ;;
        esac
    shift
    done

#   SRA to Fastq
#       -F ensures it contains only original sequence name
#       --outdir is the path to output DIR
parallel 'fastq-dump --split-files -F --gzip --outdir "${DIR}" {}' :::: ${SRA_LIST}

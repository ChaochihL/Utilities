#!/bin/bash

#PBS -l mem=3gb,nodes=1:ppn=16,walltime=24:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -u
set -o pipefail

#   Script written by Chaochih Liu
#   December 21, 2016


#   Dependencies
module load sratoolkit
module load parallel

#   Where are the .sra files located?
SRA_FILES=/panfs/roc/scratch/pmorrell/WGS_Barley

#   Where do we want our output to go?
#   As of right now, this part is hardcoded, will eventually fix
#OUT_DIR=/home/morrellp/liux1299/scratch/WBDC_Inversions_Project/WGS_raw_fastq

#   Dump fastq files into out directory
#   The -F option ensures it contains only original sequence name
#   The -O option gives it the path to output the files
#   For now, the path following -O needs to be hardcoded in order for fastq-dump to work (need to fix this)
parallel 'fastq-dump --split-files -F --gzip {} -O /home/morrellp/liux1299/scratch/WBDC_Inversions_Project/WGS_raw_fastq' ::: ${SRA_FILES}/Barke*.sra

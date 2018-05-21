#!/usr/bin/env Rscript
#   Chaochih Liu - May 21, 2018

#   This is a small script that matches accession name with geographical location, latitude, longitude, and altitude info

#   Usage: ./match_geog_loc.R [sra_accession_file] [metadata_info_file] [outfile_name.txt]

#   Where:
#   1) [sra_accession_file] is a tab-delimited file with: a) no header lines, b) 1st column containing accession names
#   2) [metadata_info_file] is a tab-delimited file with: a) header line, b) 1st column is called "Accession_Name" and contains 
#       accession names. This file contains geographical location, latitude, longitude, etc. info
#   3) [outfile_name.txt] is the full filepath to the output file including output filename and file extension .txt
#       Ex: /path/to/output_file.txt

readSraSubFiles <- function(filename) {
    df <- read.delim(
        file = filename,
        header = FALSE,
        fill = TRUE,
        na.strings = "NA"
    )
    return(df)
}

readGeogFiles <- function(filename) {
    df <- read.delim(
        file = filename,
        header = TRUE,
        fill = TRUE,
        na.strings = "NA"
    )
    return(df)
}

mergeFiles <- function(sraFile, geogFile) {
    #   Merge main file (sraFile) with file containing metadata info
    #   Merge by accession names columns
    merged <- merge(
        x = sraFile,
        y = geogFile,
        by.x = "V1",
        by.y = "Accession_Name",
        all.x = TRUE, # Rows that do not have a match will remain in dataframe
        all.y = FALSE # Rows that do not have a match with x will not be added to x
    )
    return(merged)
}

writeOutFile <- function(mergedData, outFilename) {
    write.table(
        x = mergedData,
        file = outFilename,
        sep = "\t",
        quote = FALSE,
        eol = "\n",
        col.names = TRUE,
        row.names = FALSE,
        na = ""
    )
}

main <- function() {
    #   Take in command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    sra.fp <- args[1]
    geog.fp <- args[2]
    outname <- args[3] # full filepath to output filename
    
    #   Read in files
    sra.df <- readSraSubFiles(filename = sra.fp)
    geog.df <- readGeogFiles(filename = geog.fp)
    
    #   merge files
    merged.df <- mergeFiles(sraFile = sra.df, geogFile = geog.df)
    
    #   Save to output file
    writeOutFile(mergedData = merged.df, outFilename = outname)
}

main() # Run the program

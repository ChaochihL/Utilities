#!/usr/bin/env Rscript
#   Chaochih Liu - May 21, 2018

#   This is a small script that merges SRA BioSample attributes file with metadata file

#   Usage: ./match_bioSample_attributes.R [bioSample_attributes_file] [metadata_info_file] [outfile_name.txt]

#   Where:
#   1) [bioSample_attributes_file] is a tab-delimited file with: a) no header lines, b) 1st column containing accession names
#   2) [metadata_info_file] is a tab-delimited file with: a) header line, b) 1st column is called "Accession_Name" and contains 
#       accession names. This file contains geographical location, latitude, longitude, etc. info
#   3) [outfile_name.txt] is the full filepath to the output file including output filename and file extension .txt
#       Ex: /path/to/output_file.txt

readFiles <- function(filename) {
    df <- read.delim(
        file = filename,
        header = TRUE,
        fill = TRUE,
        na.strings = "NA"
    )
    return(df)
}

mergeFiles <- function(bioSattrFile, metadataFile) {
    #   Merge main file (sraFile) with file containing metadata info
    #   Merge by accession names columns
    merged <- merge(
        x = bioSattrFile,
        y = metadataFile,
        by.x = "sample_name",
        by.y = "V1",
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
    bioSattr.fp <- args[1]
    metadata.fp <- args[2]
    outname <- args[3] # full filepath to output filename
    
    #   Read in files
    bioSattr.df <- readFiles(filename = bioSattr.fp)
    metadata.df <- readFiles(filename = metadata.fp)
    
    #   merge files
    merged.df <- mergeFiles(bioSattrFile = bioSattr.df, metadataFile = metadata.df)
    
    #   Save to output file
    writeOutFile(mergedData = merged.df, outFilename = outname)
}

main() # Run the program

#!/usr/bin/env Rscript

#   This script adds additional columns (physical positions and new Chr info for each SNP) to csv file based on matching SNP name.
#   Takes three arguments:
    #   1) VCF file
    #   2) csv file from that includes SNP name info
    #   3) Output filename
#   Script written by Chaochih Liu
#   July 21, 2016

#   To run the script: ./SNP_Position_Matching.R <positions_file> <data.csv> <output_filename>


#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   A function to read in VCF file
readVcf <- function(filename) {
    data.file <- read.table(
        file = filename, # passed as an argument
        header = FALSE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    #   For the merge we are only interested in
    #   Chr, PhysPos, and SNP_id info
    vcf.subset <- data.frame(Chr_2016 = data.file$V1,
                             PhysPos_2016 = data.file$V2,
                             SNP_id = data.file$V3)
    return(vcf.subset)
}

#   A function to read in SNP_BAC.txt file format
readCsv <- function(filename) {
    data.file <- read.csv(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    return(data.file)
}

fixNaming <- function(data.file) {
    #   Add 'X' character at the beginning of each character element
    #   regex pattern "^" is point before first character
    data.file$SNP <- sub("X", "", data.file$SNP)
    return(data.file)
}

#   A function to merge files
#   based on matching SNP names
#   If SNP names match, add position to new column in HarvEST file
mergeFile <- function(physicalData, mainFile) {
    #   Merge harvest and physical based on matches found between Query_SNP and SNP_id columns
    merged <- merge(x = mainFile,
                    y = physicalData,
                    by.x = "SNP", # Use SNP column from csv file for merge
                    by.y = "SNP_id", # USE SNP_id column from physicalData for merge
                    all.x = TRUE, # Rows that do not have a match will remain in dataframe
                    all.y = FALSE # Rows that do not have a match with x will not be added to x
                    )
    return(merged)
}

#   A function to write data to outfile
writeOutFile <- function(mergedData, outFilename) {
    write.csv(
        x = mergedData,
        file = outFilename,
        quote = FALSE,
        eol = "\n",
        col.names = TRUE,
        row.names = FALSE,
        na = "NA")
}

#   Driver function
main <- function() {
    vcfFile <- args[1] # vcf file
    csvFile <- args[2] # CSV file with column of SNP names
    outName <- args[3] # name given by user is third argument
    physical <- readVcf(filename = vcfFile) # read in physical positions
    maindata <- readCsv(filename = csvFile) # read in main csv
    maindata.no.x <- fixNaming(data.file = maindata) # remove 'X' at the beginning of SNP names
    merged <- mergeFile(physicalData = physical.no.x, mainFile = maindata.no.x) # merge physical and main csv based on matching SNP names
    writeOutFile(mergedData = merged, outFilename = outName) # write merged data to outfile
}

main() # Run the program

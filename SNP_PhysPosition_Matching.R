#!/usr/bin/env Rscript

#   This script adds additional columns to HarvEST SNP_BAC.txt with physical positions and new Chr info for each SNP.
#   Takes three arguments:
    #   1) VCF file
    #   2) SNP_BAC.txt file from HarvEST that includes genetic map/Golden Path info
    #   3) Output filename
#   Script written by Chaochih Liu
#   July 21, 2016

#   To run the script: Rscript ./SNP_Position_Matching.R <positions_file> <SNP_BAC.txt> <output_filename>


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
                             PhysPos = data.file$V2, 
                             SNP_id = data.file$V3)
    return(vcf.subset)
}

#   A function to read in SNP_BAC.txt file format
readSnpBAC <- function(filename) {
    data.file <- read.delim(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    colnames(data.file) <- c("Query_SNP", "SNP_Found", "BAC", "Seq_Source", "Node", "UCR_FPC_Contig", "Arm", "IBSC_2012_FPC_Contig", "IBSC_2012_LG", "IBSC_2012_cM", "IBSC_2012_Golden_Path", "2014_LG", "2014_cM", "2011_LG", "2011_cM", "2009_LG", "2009_cM")
    return(data.file)
}

#   A function to merge files
#   based on matching SNP names
#   If SNP names match, add position to new column in HarvEST file
mergeFile <- function(physicalData, harvestData) {
    #   Merge harvest and physical based on matches found between Query_SNP and SNP_id columns
    merged <- merge(x = harvestData,
                    y = physicalData,
                    by.x = "Query_SNP", # Use Query_SNP column from harvestData for merge
                    by.y = "SNP_id", # USE SNP_id column from physicalData for merge
                    all = FALSE # Rows that do not have a match will not be put into final dataframe
                    )
    return(merged)
}

#   A function to write data to outfile
writeOutFile <- function(mergedData, outFilename) {
    write.table(x = mergedData,
                file = outFilename,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = TRUE,
                row.names = FALSE,
                na = "NA")
}

#   Driver function
main <- function() {
    vcfFile <- args[1] # vcf file
    snpBAC <- args[2] # HarvEST Barley SNP_BAC.txt file is second arguemnt
    outName <- args[3] # name given by user is third argument
    physical <- readVcf(filename = vcfFile) # read in physical positions
    harvest <- readSnpBAC(filename = snpBAC) # read in SNP_BACs
    merged <- mergeFile(physicalData = physical, harvestData = harvest) # merge physical and harvest based on matching SNP names
    writeOutFile(mergedData = merged, outFilename = outName) # write merged data to outfile
}

main() # Run the program

#!/usr/bin/env Rscript

library(xlsx)

#   A function to read in VCF file
readXlsx <- function(filename) {
    data.file <- read.xlsx2(
        file = filename, # passed as an argument
        sheetName = "JHI-FINAL_WBDC_LIST_JAN-2017",
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    return(data.file)
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

#   A function to merge files
#   based on matching SNP names
#   If SNP names match, add position to new column in HarvEST file
mergeFile <- function(morrell_lab, jhi) {
    #   Merge harvest and physical based on matches found between Query_SNP and SNP_id columns
    merged <- merge(x = morrell_lab,
                    y = jhi,
                    by.x = "Accession_ID", # Use Query_SNP column from harvestData for merge
                    by.y = "WBDC_accession", # USE SNP_id column from physicalData for merge
                    all = FALSE # Rows that do not have a match will not be put into final dataframe
    )
    return(merged)
}

#   A function to write data to outfile
writeOutFile <- function(mergedData, outFilename) {
    write.xlsx2(x = mergedData,
              file = outFilename,
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE,
              showNA = TRUE)
}

#   Driver function
main <- function() {
    #   Take command line arguments
    #   Stores arguments into a vector
    args <- commandArgs(trailingOnly = TRUE)
    vcfFile <- args[1] # vcf file
    snpBAC <- args[2] # HarvEST Barley SNP_BAC.txt file is second arguemnt
    outName <- args[3] # name given by user is third argument
    physical <- readVcf(filename = vcfFile) # read in physical positions
    harvest <- readSnpBAC(filename = snpBAC) # read in SNP_BACs
    merged <- mergeFile(physicalData = physical, harvestData = harvest) # merge physical and harvest based on matching SNP names
    writeOutFile(mergedData = merged, outFilename = outName) # write merged data to outfile
}

main() # Run the program
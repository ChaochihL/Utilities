#!/usr/bin/env Rscript

#   This script removes all rows with empty cells in Query_SNP column and missing info for PhysPos column.
#   Takes one argument:
    #   1) SNP_BAC.txt file that includes PhysPos and Chr_2016 columns
#   Script written by Chaochih Liu
#   September 20, 2016

#   To run the script: Rscript ./snpBAC_NA_filtering.R <SNP_BAC.txt>


#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   A function to read in SNP_BAC.txt file format that includes
#   physical positions and new chromosome info
readFileSnpBAC <- function(filename) {
    data.file <- read.delim(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    return(data.file)
}

#   A function to remove rows with empty cells or missing PhysPos info
removeEmptyAndNA <- function(data.file) {
    tmpEmptyRemoved <- data.file[!(data.file$Query_SNP == ""), ] # if Query_SNP column contains empty cell, remove row
    filtered <- tmpEmptyRemoved[!is.na(tmpEmptyRemoved$PhysPos), ] # if PhysPos column contains NA value, remove row
    return(filtered)
}

#   Write file to outfile
writeOutFile <- function(filename, filtered) {
    inputName <- unlist(strsplit(x = filename, split = ".txt"))
    outputName <- paste(inputName, "filtered.txt", sep = "_")
    write.table(x = filtered,
                file = outputName,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = TRUE,
                row.names = FALSE)
}

#   Driver function
main <- function() {
    snpBACPhysPos <- args[1] # SNP_BAC.txt file with Query_SNP and PhysPos columns
    dataToFilter <- readFileSnpBAC(filename = snpBACPhysPos) # read in SNP_BAC.txt file
    finalFiltered <- removeEmptyAndNA(data.file = dataToFilter) # data that has rows with empty cells or NA values for PhysPos removed
    writeOutFile(filename = snpBACPhysPos, filtered = finalFiltered) # write filtered data to outfile
}

main() # Run the program
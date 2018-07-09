#!/usr/bin/env Rscript

#   This script converts columns in dataframe to genotype format accepted by LDheatmap R package.
#   Depends on R package: genetics
#   Script written by Chaochih Liu
#   November 28, 2016

library(genetics)

#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   Read in genotype data
#   Input data frame should be sorted by SNP names
readGenoFile <- function(filename) {
    #   we want the SNPs as columns & sample names (i.e. WBDC) as rows
    #   the file we read in has sample names (i.e. WBDC) as columns & SNPs as rows
    genoData <- t(read.delim(file = filename, 
                             row.names = 1, # include row names
                             header = TRUE)) # include column names
    return(genoData)
}

makeGeno <- function(genoData) {
    #   Convert genotype format to format that works with LDheatmap function
    #   Example: flrom AA to A/A
    tmp.geno <- makeGenotypes(data = genoData, sep = "")
    return(tmp.geno)
}

#   Write out file
writeOutFile <- function(filename, gData) {
    inputName <- unlist(strsplit(x = filename, split = ".txt"))
    outputName <- paste0(inputName, "_MG.txt")
    write.table(x = gData,
                file = outputName,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = TRUE,
                row.names = TRUE)
}

#   Driver function
main <- function() {
    genoFile <- args[1] # genotype data
    geno <- readGenoFile(filename = genoFile)
    tmp.geno <- makeGeno(genoData = geno)
    writeOutFile(filename = genoFile, gData = tmp.geno)
}

main() # Run the program

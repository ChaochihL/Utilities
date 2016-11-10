#!/usr/bin/env Rscript

#   This script converts our HarvEST SNP_BAC.txt file with PhysPos/Chr_2016 info into BED file
#   to send to bedops.
#   NOTE: this is not a true BED file, it includes an extra column with SNP ID

#   Takes one argument:
    #   1) SNP_BAC.txt file that includes Query_SNP, PhysPos, and Chr_2016 columns
#   Script written by Chaochih Liu
#   November 10, 2016

args <- commandArgs(trailingOnly = TRUE)

#   Read in file containing SNP names, physical position, and chromosome info
readSnpBac <- function(filename) {
    tmp.data.file <- read.delim(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    data.file <- data.frame(Chr_2016 = tmp.data.file$Chr_2016, 
                            PhysPos = tmp.data.file$PhysPos, 
                            Query_SNP = tmp.data.file$Query_SNP)
    return(data.file)
}

#   Sort, filter, and generate BED file
processData <- function(data.file) {
    tmp.uniq <- unique(x = data.file) # Remove duplicate PhysPos/SNPs
    tmp.uniq.ordered <- tmp.uniq[order(tmp.uniq$PhysPos), ] # sort by PhysPos
    #   Subtract 1 from Physical position & add as new column
    #   BED files are 0-indexed
    tmp.uniq.ordered["Start"] <- (tmp.uniq.ordered$PhysPos - 1)
    #   Reorder output dataframe
    final.df <- data.frame(chr = tmp.uniq.ordered$Chr_2016,
                           start = tmp.uniq.ordered$Start,
                           end = tmp.uniq.ordered$PhysPos,
                           ID = tmp.uniq.ordered$Query_SNP)
    return(final.df)
}

#   Save new dataframe as BED file
#   except with extra column for SNP ID at the end
writeOutBed <- function(filename, final.df) {
    inputName <- unlist(strsplit(x = filename, split = ".txt"))
    outputName <- paste(inputName, ".bed", sep = "")
    write.table(x = final.df,
                file = outputName,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = FALSE,
                row.names = FALSE)
}

#   Driver function
main <- function() {
    snpBacFile <- args[1] # SNP_BAC.txt file with Query_SNP, PhysPos, and Chr_2016 columns
    snpBac.df <- readSnpBac(filename = snpBacFile) # read in data
    bedFile <- processData(data.file = snpBac.df) # create BED file
    writeOutBed(filename = snpBacFile, final.df = bedFile) # save data to file
}

main() # Run the program

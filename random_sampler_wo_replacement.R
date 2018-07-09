#!/usr/bin/env Rscript

#   This script randomly selects SNP names from a given list
#   Script written by Chaochih Liu
#   June 17, 2016

#   To run the script: ./random_sampler.R [number_to_sample] [sample_list.txt]

#   User provided input arguments are:
#       1) number_to_sample - how many do we want to sample?
#       2) sample_list.txt is a file with a single column of values and NO header line

#   Read in sample list and create data frame
readFile <- function(filename) {
    sampleList <- read.delim(
        file = filename,
        header = FALSE # file does not have header, don't want 1st row to turn into header
    )
    return(sampleList)
}

#   Random sampler
randomSampler <- function(sampleList, numberOfDraws) {
    set.seed(Sys.time())
    randomDraws <- sort(
        sample(
            x = as.vector(sampleList[, 1]), # sort samples from lowest to highest number 
            size = numberOfDraws, # Number of samples we want to draw
            replace = FALSE
        )
    )
    return(randomDraws)
}

#   Write file to outfile
writeOutFile <- function(randomDraws, filename, numberOfDraws) {
    inputName <- unlist(strsplit(x = filename, split = ".txt"))
    outputName <- paste(inputName, as.character(x = numberOfDraws), "randomSamples.txt", sep = "_")
    write.table(
        x = randomDraws,
        file = outputName,
        quote = FALSE,
        sep = "\t",
        eol = "\n",
        col.names = FALSE,
        row.names = FALSE
    )
}

#   Run the script
runScript <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    nSelect <- as.numeric(args[1]) # Number of samples to draw
    inputData <- args[2] # filename
    sampleNames <- readFile(inputData)
    samplesSelected <- randomSampler(sampleNames, nSelect)
    writeOutFile(randomDraws = samplesSelected, filename = inputData, numberOfDraws = nSelect)
}

runScript() # Run program

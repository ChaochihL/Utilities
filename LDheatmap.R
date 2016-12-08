#!/usr/bin/env Rscript

#   This is a script that generates LD heatmap
#   The concept of this script was drawn from Peter Morrell's script: 
    #   https://github.com/pmorrell/Utilities/blob/master/LDheatmap/CAP_heatmap_7H.r
#   Written by Chaochih Liu
#   September 26, 2016

#   Required arguments:
    #   1) genoData.txt: genotype data frame that has SNPs sorted
    #   2) physPos.txt: file that includes two columns - "Query_SNP" and "PhysPos"
    #   3) heatmap plot name: this will be the title of the heatmap
    #   4) outName: outfile prefix (do not include file extension)

#   To run: Rscript ./LDheatmap.R <genoData.txt> <physPos.txt> <Plot Name> <Out File Name>


library(LDheatmap)
library(genetics)
library(RColorBrewer)
library(grDevices)
library(chopsticks)
require(chopsticks) # used for LDheatmap function

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

#   Read in file containing SNP names and physical position info
#   Input file should already be sorted by Query_SNP
readPhysPos <- function(filename) {
    tmp.data.file <- read.delim(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    data.file <- data.frame(Query_SNP = tmp.data.file$Query_SNP, PhysPos = tmp.data.file$PhysPos)
    return(data.file)
}

#   Fix naming of markers to match genoData marker naming
fixNaming <- function(data.file) {
    #   Add 'X' character at the beginning of each character element
    #   regex pattern "^" is point before first character
    data.file$Query_SNP <- sub("^", "X", data.file$Query_SNP)
    uniq.file <- unique(x = data.file) # keep unique SNP names only
    return(uniq.file)
}

#   Create names used for heatmap plot
hm.plotNames <- function(data.file) {
    #   Keep unique SNP names only
    hm.naming <- unique(x = data.file)
    return(hm.naming)
}

#   Generate heatmap for r^2 values
hm.r2 <- function(genoData, PhysPos, plotName, outName) {
    #   Convert genotype format AA to A/A
    #   This format works with LDheatmap function
    tmp.genos <- makeGenotypes(data = genoData, sep = "")
    #   Generate heatmap color palette
    heatmapColors <- brewer.pal(n = 9, name = "YlOrRd")
    #   Output file naming
    outputName <- paste("HM", "r2", sep = "_")
    name <- paste0(outName, "_", outputName, ".png")
    #png(filename = name)
    png(filename = name)
    #   Run LDheatmap
    heatmap.r2 <- LDheatmap(gdat = tmp.genos,
                            genetic.distances = PhysPos$PhysPos,
                            distances = "physical",
                            LDmeasure = "r",
                            add.map = FALSE,
                            add.key = TRUE,
                            SNP.name = PhysPos$Query_SNP,
                            flip = FALSE,
                            color = heatmapColors,
                            title = plotName)
    dev.off()
    return(heatmap.r2)
}

#   Generate heatmap for D' values
hm.Dprime <- function(genoData, PhysPos, plotName, outName) {
    #   Convert genotype format AA to A/A
    #   This format works with LDheatmap function
    tmp.genos <- makeGenotypes(data = genoData, sep = "")
    #   Generate heatmap color palette
    heatmapColors <- brewer.pal(n = 9, name = "YlOrRd")
    #   Output file naming
    outputName <- paste("HM", "Dprime", sep = "_")
    name <- paste0(outName, "_", outputName, ".png")
    png(filename = name)
    #   Generate D' plot
    heatmap.D <- LDheatmap(gdat = tmp.genos,
                           genetic.distances = PhysPos$PhysPos,
                           distances = "physical",
                           LDmeasure = "D'",
                           add.map = FALSE,
                           add.key = TRUE,
                           flip = FALSE,
                           color = heatmapColors,
                           title = plotName)
    dev.off()
    return(heatmap.D)
}

#   Save heatmap to .png file
outFile <- function(outName, heatmap, LDcalc) {
    #inputName <- unlist(strsplit(x = as.character(quote(LDcalc)), split = "heatmap."))[2]
    outputName <- paste("HM", LDcalc, sep = "_")
    #   Name output .png file
    name <- paste0(outName, "_", outputName, ".txt")
    heatmapDF <- as.data.frame(x = heatmap$LDmatrix, row.names = NULL)
    write.table(x = heatmapDF,
                file = name,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = NA, # To prevent top row from shifting due to empty first cell
                row.names = TRUE) # Want to keep SNP names
}

#   Driver function
main <- function() {
    genoData <- args[1] # 1) genotype data frame that has SNPs sorted
    physPosData <- args[2] # 2) file that includes two columns - "Query_SNP" and "PhysPos"
    plotName <- args[3] # 3) heatmap plot name
    outName <- args [4] # 4) outfile prefix (do not include file extension)
    genoFile <- readGenoFile(filename = genoData)
    physPos <- readPhysPos(filename = physPosData)
    X.names <- fixNaming(data.file = physPos) # marker names with 'X' as first character
    snpNames <- hm.plotNames(data.file = physPos) # marker names to use for heatmap plot
    #   r^2 heatmap
    plot.r2 <- hm.r2(genoData = genoFile, 
                     PhysPos = X.names, 
                     plotName = plotName,
                     outName = outName)
    #   D' heatmap
    plot.D <- hm.Dprime(genoData = genoFile,
                        PhysPos = X.names,
                        plotName = plotName,
                        outName = outName)
    #   Save r2 plot to .png file
    outFile(outName = outName,
            heatmap = plot.r2,
            LDcalc = "r2")
    #   Save D' plot to .png file
    outFile(outName = outName,
            heatmap = plot.D,
            LDcalc = "D_prime")
}

main() # Run the program

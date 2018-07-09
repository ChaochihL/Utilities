#!/usr/bin/env Rscript

#   This is a script that generates LD heatmap: r2 and D'
#   The concept of this script was drawn from Peter Morrell's script:
    #   https://github.com/pmorrell/Utilities/blob/master/LDheatmap/CAP_heatmap_7H.r
#   Written by Chaochih Liu
#   September 26, 2016

#   Required arguments:
    #   1) genoData.txt: genotype data frame that has SNPs sorted by SNP names
    #   2) physPos.txt: file that includes two columns named - "Query_SNP" and "PhysPos"
    #   3) heatmap plot name: this will be the title of the heatmap figure
    #   4) outPrefix: outfile prefix (do not include file extension)
    #   5) outDir: full path to where output files go
    #   6) includeSnpName: only takes in two args - include/exclude
    #   7) # of individuals: total # of individuals (i.e. total # of WBDC samples)
    #   8) missing data threshold: proportion of missing data to use for filtering
    #       i.e. if we want to filter out greater than 20% missing data,
    #       our argument input would be 0.20

#   To run: ./LDheatmap.R [genoData.txt] [physPos.txt] [Plot Name] [Out File Prefix] [out directory] [include/exclude] [# of individuals] [missing data threshold]

library(LDheatmap)
library(genetics)
library(RColorBrewer)
library(grDevices)
require(chopsticks) # used for LDheatmap function

#   Read in genotype data with input data frame sorted by SNP names
readGenoFile <- function(filename) {
    #   we want the SNPs as columns & sample names (i.e. WBDC) as rows
    #   the file we read in has sample names (i.e. WBDC) as columns & SNPs as rows
    genoData <- t(read.delim(file = filename,
                             row.names = 1, # include row names
                             header = TRUE,
                             na.strings = "NN") # 'NN' is missing data, convert to NA
    )
    return(genoData)
}

#   Read in file containing SNP names and physical position info
#   Input file should already be sorted by Query_SNP
readPhysPos <- function(filename) {
    tmp.data.file <- read.delim(
        file = filename,
        header = TRUE,
        fill = TRUE,
        na.strings = "NA"
    )
    data.file <- data.frame(Query_SNP = tmp.data.file$Query_SNP,
                            PhysPos = tmp.data.file$PhysPos)
    return(data.file)
}

#   Convert genotype format to LDheatmap function compatible format
#   Example: from AA to A/A
makeGeno <- function(genoData) {
    cat("Running makeGenotypes...", sep = "\n")
    #   This also converts NA to N/A
    compatible.df <- makeGenotypes(data = genoData, sep = "")
    cat("Done running makeGenotypes...", sep = "\n")
    return(compatible.df)
}

#   Remove empty columns
findAllNA <- function(df.column, n.individuals) {
    #   Test if "N/A" present in columns
    #   If number of "N/A" is equal to n.individuals (meaning all empty cells),
    #   remove marker
    search.all.na <- sum(grepl(pattern = "N/A", x = df.column))
    finding <- search.all.na == n.individuals
    return(finding)
}

#   Find incompatible markers in genotype data
#   Paul Hoffman assisted in writing code for this function
findIncompatible <- function(df.column, n.individuals) {
    #   Use grepl to return logical vector of TRUE/FALSE for being compatible
    #   If compatible (i.e. A/A) or N/A, will return TRUE
    #   If incompatible (i.e. AA), will return FALSE
    searches <- grepl(pattern = "[ACTGN]/[ACTGN]", x = df.column)
    #   Sum all the N/A or incompatible data
    #   Invert all TRUEs and FALSEs to better find failures using sum()
    #   Now, all TRUEs are incompatible
    #   And all FALSEs are compatible
    search.summary <- sum(!searches)
    #   This works because all trues inverted become zero
    #   So, a single fail means that we get a sum of greater than zero
    search.bool <- search.summary == n.individuals
    return(search.bool)
}

#   Filter out greater than n% missing data (specified in input argument)
filterMissing <- function(df.column, n.missing) {
    #   If number of NA's is greater than n.missing (our threshold),
    #   we will filter out those markers
    #   If "N/A" pattern is found, it will return TRUE
    na.searches <- grepl(pattern = "N/A", x = df.column)
    #   Are number of NAs greater than threshold?
    #   If FALSE, keep the marker
    #   If TRUE, filter out the marker
    na.summary <- sum(na.searches) | sum(is.na(df.column)) > n.missing
    return(na.summary)
}

#   Generate heatmap for r^2 values
hm.r2 <- function(genoData, PhysPos, plotName, snpName, outName, directory) {
    heatmapColors <- rev(brewer.pal(n = 9, name = "YlOrRd"))
    outputName <- paste("HM", "r2", sep = "-")
    pdf(file = paste0(directory, "/", outName, "-", outputName, ".pdf"))
    heatmap.r2 <- LDheatmap(gdat = genoData,
                            genetic.distances = PhysPos$PhysPos,
                            distances = "physical",
                            LDmeasure = "r",
                            add.map = FALSE,
                            add.key = TRUE,
                            SNP.name = snpName,
                            flip = FALSE,
                            color = heatmapColors,
                            title = plotName)
    dev.off()
    return(heatmap.r2)
}

#   Generate heatmap for D' values
hm.Dprime <- function(genoData, PhysPos, plotName, snpName, outName, directory) {
    heatmapColors <- rev(brewer.pal(n = 9, name = "YlOrRd"))
    outputName <- paste("HM", "Dprime", sep = "-")
    pdf(file = paste0(directory, "/", outName, "-", outputName, ".pdf"))
    heatmap.D <- LDheatmap(gdat = genoData,
                           genetic.distances = PhysPos$PhysPos,
                           distances = "physical",
                           LDmeasure = "D'",
                           add.map = FALSE,
                           add.key = TRUE,
                           SNP.name = snpName,
                           flip = FALSE,
                           color = heatmapColors,
                           title = plotName)
    dev.off()
    return(heatmap.D)
}

#   Save info to spreadsheet
outCsv <- function(df, rowNames, outName) {
    write.csv(x = df,
              file = outName,
              row.names = rowNames)
}

#   Write data to output file
outCompatibleFile <- function(df, outName, outDirectory) {
    name <- paste0(outDirectory, "/", outName, ".txt")
    write.table(x = df,
                file = name,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = TRUE, # To prevent top row from shifting due to empty first cell
                row.names = TRUE) # Want to keep SNP names
}

#   Save heatmap to .pdf file
outFile <- function(outName, directory, heatmap, LDcalc) {
    outputName <- paste("HM", LDcalc, sep = "_")
    name <- paste0(directory, "/", outName, "_", outputName, ".txt")
    heatmap.df <- as.data.frame(x = heatmap$LDmatrix, row.names = NULL)
    write.table(x = heatmap.df,
                file = name,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = NA, # To prevent top row from shifting due to empty first cell
                row.names = TRUE) # Want to keep SNP names
}

#   Driver function
main <- function() {
    #   Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    genoData <- args[1] # 1) genotype data frame that has SNPs sorted
    physPosData <- args[2] # 2) file that includes two columns - "Query_SNP" and "PhysPos"
    plotName <- args[3] # 3) heatmap plot name
    outPrefix <- args [4] # 4) outfile prefix (do not include file extension)
    outDir <- args[5] # 5) where do our output files go?
    includeSnpName <- args[6] # 6) Do we want to include SNP names in our plot? (include/exclude)
    #   Filter out n% missing data
    n.individuals <- as.numeric(args[7]) # Our dataset has 284 individuals (i.e. WBDC)
    p.missing <- as.numeric(args[8]) # What percent missing data threshold are we using?
    indv.missing <- round(n.individuals * p.missing) # The number of missing individuals we will use as threshold

    #   Read in genotype data and SNP_BAC.txt file
    genoFile <- readGenoFile(filename = genoData)
    physPos <- readPhysPos(filename = physPosData)
    X.names <- physPos # We don't need to fix names in this dataset, just use physPos

    #   Create output and plot names
    outname <- paste0(outDir, "/", outPrefix, "_SNP_info")

    #   Convert data to compatible format
    geno.converted <- makeGeno(genoData = genoFile)

    cat("Removing empty columns filled with NA's...", sep = "\n")
    results.na <- apply(X = geno.converted,
                        MARGIN = 2,
                        FUN = findAllNA,
                        n.individuals = n.individuals)
    #   Keep a record of empty columns
    no.data.cols <- names(x = which(x = results.na))
    cat("Saving samples with empty columns to spreadsheet.", sep = "\n")
    tryCatch({
        outCsv(df = no.data.cols,
                rowNames = FALSE,
                outName = paste0(outname, "-empty_cols.csv"))
    }, error = function(e) {
        cat("No empty columns found.", sep = "\n")
        }
    )
    #   Data frame excluding emtpy columns
    results.na.removed <- geno.converted[, !results.na]

    #   Filter out incompatible genotype columns
    cat("Removing incompatible columns...", sep = "\n")
    results <- apply(X = results.na.removed,
                     MARGIN = 2, # Applied over columns
                     FUN = findIncompatible,
                     n.individuals = n.individuals)
    #   Keep a record of column names that failed
    failed.samples <- names(x = which(x = results))
    cat("Saving failed samples to spreadsheet.", sep = "\n")
    tryCatch({
        outCsv(df = failed.samples,
                rowNames = FALSE,
                outName = paste0(outname, "-failed_snps.csv"))
    }, error = function(e) {
        cat("No failed samples.", sep = "\n")
        }
    )
    #   Remove markers with incompatible genotype columns
    incomp.cols.removed <- results.na.removed[, !results]

    #   Filter out data with greater than n% missing data
    cat("Filtering out data with greater than:", sep = "\n")
    cat(p.missing * 100, "% missing data")
    results.filt <- apply(X = incomp.cols.removed,
                          MARGIN = 2, # Applied over columns
                          FUN = filterMissing,
                          n.missing = indv.missing)
    #   Keep a record of markers that were filtered out
    miss.data.cols <- names(x = which(x = results.filt))

    cat("Saving filtered out markers to spreadsheet.", sep = "\n")
    tryCatch({
        #   Markers that had greater than n% missing data will get saved
        outCsv(df = miss.data.cols,
               rowNames = FALSE,
               outName = paste0(outname, "-missing_data_cols.csv"))
    }, error = function(e) {
        cat("No markers with greater than:", sep = "\n")
        cat(p.missing *100, "% missing data ")
        cat("found.", sep = "\n")
        }
    )

    #   Keep compatible data after removing markers with more than n% missing data
    pass.samples <- incomp.cols.removed[, !results.filt]
    cat("Saving compatible data used for heatmap to file.", sep = "\n")
    outCompatibleFile(
        df = pass.samples,
        outName = paste0(outPrefix, "_compatibleSnps"),
        outDirectory = outDir
    )

    #   Remove SNPs with empty columns, missing data > threshold, and failed SNPs from PhysPos file
    X.names.noEmpty <- X.names[!(X.names$Query_SNP %in% no.data.cols), ]
    X.names.failed <- X.names.noEmpty[!(X.names.noEmpty$Query_SNP %in% failed.samples), ]
    X.names.filtered <- X.names.failed[!(X.names.failed$Query_SNP %in% miss.data.cols), ]

    #   Do we include SNP names or exclude SNP names in our heatmap plot?
    if(includeSnpName == "exclude") {
        cat("Exclude SNP names from plots.", sep = "\n")
        snpname <- FALSE
    } else {
        cat("Include SNP names in plots.", sep = "\n")
        snpname <- X.names.filtered$Query_SNP
    }

    cat("LDheatmap - r2 starting...", sep = "\n")
    plot.r2 <- hm.r2(genoData = pass.samples,
                     PhysPos = X.names.filtered,
                     plotName = plotName,
                     snpName = snpname,
                     outName = outPrefix,
                     directory = outDir)
    cat("LDheatmap - r2 done.", sep = "\n")

    cat("LDheatmap - D' starting...", sep = "\n")
    plot.D <- hm.Dprime(genoData = pass.samples,
                        PhysPos = X.names.filtered,
                        plotName = plotName,
                        snpName = snpname,
                        outName = outPrefix,
                        directory = outDir)
    cat("LDheatmap - D' done", sep = "\n")

    cat("Saving files to out directory...", sep = "\n")
    outFile(outName = outPrefix,
            directory = outDir,
            heatmap = plot.r2,
            LDcalc = "r2")
    outFile(outName = outPrefix,
            directory = outDir,
            heatmap = plot.D,
            LDcalc = "Dprime")
    cat("Done.", sep = "\n")
}

#   Run the program
main()

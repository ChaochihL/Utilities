#!/usr/bin/env Rscript

# This script takes output .hist.all.txt files from sequence_handling coverage_mapping and
# generates coverage summary plots.

# Usage:
#   ./coverage_plotting.R [hist_dir] [project] [out_dir]

# Where:
#   1) [hist_dir] is the full filepath to directory containing .hist.all.txt files
#       Note: Exclude the last '/' at the end of the filepath. Script adds the slash at the end.
#   2) [project] is the project name, this will be used as prefix for plot filenames.
#   3) [out_dir] is the full directory to where our plots should output to.
#       Note: Exclude the last '/' at the end of the filepath. Script adds the slash at the end.

# IMPORTANT: This script assumes you have run the following command:
#   to pull out histogram summarizing coverage among "all" features in A
#   grep ^all ${out_dir}/Histograms/${sampleName}.hist > ${out_dir}/Histograms/${sampleName}.hist.all.txt

# Code adapted from:
# https://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html

library(RColorBrewer)
library(ggplot2)
library(data.table)

readHistFiles <- function(hist_dir, fp) {
    # Read in coverage histogram file
    cov <- list()
    # Using fread to read in file instead of read.table or read.delim2
    # because fread is significantly faster than both
    for (i in 1:length(fp)) {
        cov[[i]] <- fread(file = fp[i], header = FALSE, sep = "\t")
    }
    return(cov)
}

plotCoverage <- function(cov, cov_cumulative, max_depth, cols_more, labs, main_title) {
    # Create plot area, but do not plot anything. Add gridlines and axis labels.
    plot(
        as.numeric(as.character(unlist(cov[[1]][2:max_depth, 2]))),
        as.numeric(as.character(unlist(cov_cumulative[[1]][1:max_depth-1]))),
        type='n', xlab="Depth", ylab=expression("Fraction of bases" >= "depth"),
        ylim=c(0,1.0), main=main_title
    )
    abline(v = 20, col = "gray60")
    abline(v = 50, col = "gray60")
    abline(v = 80, col = "gray60")
    abline(v = 100, col = "gray60")
    abline(h = 0.50, col = "gray60")
    abline(h = 0.90, col = "gray60")
    axis(1, at=c(20,50,80), labels=c(20,50,80))
    axis(2, at=c(0.90), labels=c(0.90))
    axis(2, at=c(0.50), labels=c(0.50))

    # Actually plot the data for each alignment (stored in lists)
    for (i in 1:length(cov)) {
        points(
            as.numeric(as.character(unlist(cov[[i]][2:max_depth, 2]))),
            as.numeric(as.character(unlist(cov_cumulative[[i]][1:max_depth-1]))),
            type = 'l',
            lwd = 2.5,
            col = cols_more[[i]]
        )
    }

    # Add a legend using the nice sample labeles rather than the full filenames.
    legend("topright", legend=labs, col=cols_more, lty=1, lwd=4, ncol = 4, cex = 0.6)
}

main <- function() {
    # Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    # User provided command line arguments
    hist_dir <- args[1]
    project <- args[2]
    out_dir <- args[3]

    # Prepping filepaths
    files <- list.files(path = hist_dir, pattern = ".hist")
    fp <- paste0(hist_dir, "/", files)

    # Read in files and do some processing
    cov <- readHistFiles(hist_dir, fp)

    # Get cumulative coverage for each alignment
    cov_cumulative <- list()
    for (i in 1:length(fp)) {
        cov_cumulative[[i]] <- 1-cumsum(cov[[i]][,5])
    }

    # Prep legend labels for each sample
    sample_names <- basename(fp)
    labs <- gsub(pattern = ".hist", replacement = "", x = sample_names)

    # Pick some colors
    # "Dark2" color palette has a limit of 8 colors
    cols <- brewer.pal(8, "Dark2")
    # Here is a trick to add more than 8 colors
    cols_more <- colorRampPalette(cols)(length(cov))

    # Save plot to PDF
    pdf(file = paste0(out_dir, "/", project, "_coverage_plot-maxdepth.pdf"),
        width = 10, height = 7)
    # Plot all the way until maximum depth
    plotCoverage(cov, cov_cumulative, max_depth = length(cov[[1]]$V2), cols_more, labs,
                 main_title = "Coverage")
    dev.off()

    # Save plot to PDF
    pdf(file = paste0(out_dir, "/", project, "_coverage_plot-depth120.pdf"),
        width = 10, height = 7)
    # Plot only until depth of 120
    plotCoverage(cov, cov_cumulative, max_depth = 120, cols_more, labs,
                 main_title = "Coverage")
    dev.off()
}

main() # Run the program

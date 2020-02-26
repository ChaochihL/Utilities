#!/usr/bin/env Rscript

# This script randomly selects SNP names from a given list
# Currently supports reading in gzipped file
# Chaochih Liu - June 17, 2016

# To run the script: ./random_sampler.R [number_to_sample] [sample_list.txt] [out_dir]

# User provided input arguments are:
#   [sample_list.txt] is a file with a single column of values and NO header line
#   [number_to_sample] is how many do we want to sample?
#   [out_dir] is the full path to where we want to output our sampling

read_file <- function(filename) {
    zfile <- gzfile(filename, 'rt')
    sample_list <- read.delim(
        file = zfile,
        header = TRUE
    )
    return(sample_list)
}

# Randomly sample rows from data frame
random_sampler <- function(dat, num_draws) {
    set.seed(Sys.time())
    random_draws <- dat[sample(x = nrow(dat), size = num_draws, replace = FALSE), ]
    sorted_draws <- random_draws[order(random_draws$marker), ]
    return(sorted_draws)
}

write_out_file <- function(random_draws, filename, num_draws, out_dir) {
    input_name <- unlist(strsplit(basename(filename), split = ".beagle.gz"))
    output_name <- paste(out_dir,"/",input_name,"_",as.character(x = num_draws), "_randomSamples.txt", sep = "")
    write.table(
        x = random_draws,
        file = output_name,
        quote = FALSE,
        sep = "\t",
        eol = "\n",
        col.names = TRUE,
        row.names = FALSE
    )
}

main <- function() {
    # Driver function
    args <- commandArgs(trailingOnly = TRUE)
    input_data <- args[1] # filename
    n_select <- as.numeric(args[2]) # Number of samples to draw
    output_dir <- args[3]
    sample_names <- read_file(input_data)
    samples_selected <- random_sampler(sample_names, n_select)
    write_out_file(random_draws = samples_selected, filename = input_data, num_draws = n_select, out_dir = output_dir)
}

main() # Run program

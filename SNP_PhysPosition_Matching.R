#   This script adds an additional physical positions column to a given file. Physical positions are added based on matching SNP names between two files.
#   Takes three arguments:
#   1) positions_file.txt: includes SNP names and Physical positions columns
#   2) snp_info.txt: new column of physical positions will be added to this file (i.e. HarvEST SNP_BAC.txt)
#   3) Output filename
#   Script written by Chaochih Liu
#   July 21, 2016

#   To run the script: Rscript ./SNP_PhysPosition_Matching.R <positions_file.txt> <snp_info.txt> <output_filename>

#   NOTE: all column names must use underscore NOT a space!


#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   A function to read in positions_file.txt format
readFilePos <- function(filename) {
    data.file <- read.table(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    return(data.file)
}

#   A function to read in snp_info.txt file format
readFileSnpInfo <- function(filename) {
    data.file <- read.delim(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    return(data.file)
}

#   A function to merge files
#   based on matching SNP names
#   If SNP names match, add position to new column in snp_info.txt file
mergeFile <- function(physicalData, snpInfoData) {
    #   Merge harvest and physical based on matches found between Query_SNP and SNP_id columns
    merged <- merge(x = snpInfoData,
                    y = physicalData,
                    #by.x = "Query_SNP", # Use SNP_name column from snp_info.txt for merge
                    by.x = "SNP_Name",
                    by.y = "SNP_id", # USE SNP_id column from physicalData for merge
                    all = FALSE # Rows that do not have a match will not be put into final dataframe
    )
    return(merged)
}

#   A function to write data to outfile
writeOutFile <- function(mergedData, outFilename) {
    #colnames(mergedData) <- c("Query_SNP", "SNP_Found", "BAC", "Seq_Source", "Node", "UCR_FPC_Contig", "Arm", "IBSC_2012_FPC_Contig", "IBSC_2012_LG", "IBSC_2012_cM", "IBSC_2012_Golden_Path", "2014_LG", "2014_cM", "2011_LG", "2011_cM", "2009_LG", "2009_cM", "Phys_Position")
    #colnames(mergedData) <- c("SNP_Name", "Chr", "cM", "FST", "GenBank_ID", "Gene_Short_Name", "Position", "Silent", "Physical Position")
    #colnames(mergedData) <- c("SNP_name", "Start", "End")
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
    physPositionFile <- args[1] # file with physical positions is first argument
    snpInfoFile <- args[2] # snp_info.txt file is second argument
    outName <- args[3] # name given by user is third argument
    physical <- readFilePos(filename = physPositionFile) # read in physical positions
    snpInfo <- readFileSnpInfo(filename = snpInfoFile) # read in snp_info.txt
    merged <- mergeFile(physicalData = physical, snpInfoData = snpInfo) # merge physical and snpInfo based on matching SNP names
    writeOutFile(mergedData = merged, outFilename = outName) # write merged data to outfile
}

main() # Run the program

#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 3) {
  stop("Usage: Rscript run_bayespeak.R <signal_file> <control_file> <result_file>")
}

# Extract arguments
signal <- args[1]
control <- args[2]
result <- args[3]

# Load required library
library(BayesPeak)

# Run bayespeak analysis
raw.output <- bayespeak(signal, control)

# Summarize peaks
output <- summarize.peaks(raw.output, method = "lowerbound")

# Load GenomicRanges
library(GenomicRanges)

# Convert RangedData to GRanges
gr_output <- as(output, "GRanges")

# Convert GRanges to data frame
df_output <- as.data.frame(gr_output)

# Write to CSV
write.csv(df_output, file = result, quote = FALSE, row.names = FALSE)

# Print confirmation message
cat("Analysis complete. Results written to:", result, "\n")
#!/usr/bin/env Rscript

# Run Subread to map reads to the genome.
#
# To submit on RCC Midway:
#
#   sbatch --mem=12g --partition=broadwl run-subread.R <nthreads> <index> <readfile1> <readfile2> <output_file>

# Input ------------------------------------------------------------------------

# Obtain name of fastq file passed at the command line
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
nthreads <- as.integer(args[1])
index <- args[2]
readfile1 <- args[3]
readfile2 <- args[4]
output_file <- args[5]
stopifnot(file.exists(readfile1, readfile2))

# Setup ------------------------------------------------------------------------

library("Rsubread")

# Create output directory if needed
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Map reads with Subread -------------------------------------------------------

align(index = index, readfile1 = readfile1, readfile2 = readfile2,
      output_file = output_file, nthreads = nthreads)

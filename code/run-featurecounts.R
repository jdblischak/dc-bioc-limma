#!/usr/bin/env Rscript

# Run featureCounts to counts reads per gene.
#
# To submit on RCC Midway:
#
#   sbatch --mem=4G --nodes=1 --tasks-per-node=4 --partition=broadwl run-featurecounts.R <nthreads> <exons> <file.bam> <output_file>

# Input ------------------------------------------------------------------------

# Obtain name of BAM file passed at the command line
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
nthreads <- as.integer(args[1])
# Name of annotation file in SAF format
saf <- args[2]
stopifnot(file.exists(saf))
bam <- args[3]
stopifnot(file.exists(bam))
output_file <- args[4]

# Setup ------------------------------------------------------------------------

library("Rsubread")

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Counts reads with featureCounts ----------------------------------------------

# Require that both pairs of the fragment map to same chromosome to be counted
counts <- featureCounts(files = bam, annot.ext = saf, isPairedEnd = TRUE,
                        requireBothEndsMapped = TRUE,
                        countChimericFragments = FALSE, nthreads = nthreads)

counts_out <- data.frame(rownames(counts$counts),
                         counts$counts[, 1])
colnames(counts_out) <- c("gene", basename(bam))

write.table(counts_out, file = output_file, quote = FALSE, sep = "\t",
            row.names = FALSE)

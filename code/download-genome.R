#!/usr/bin/env Rscript

# Download and index mouse genome for mapping with Subread.
#
# To submit on RCC Midway:
#
#   sbatch --mem=12G --partition=broadwl download-genome.R outdir
#
# Notes:
#
#  * The indexing step requires at least 8 GB of memory.
#
#  * The available releases on the Ensembl FTP site can be viewed at
#    ftp://ftp.ensembl.org/pub/

# Input ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
# Directory to save genome fasta files
# outdir <- "/project2/gilad/jdblischak/dc/genome"
outdir <- args[1]

# Chromsomes for mapping
chroms <- c(1:19, "X", "Y", "MT")

# Ensembl release to use for annotation
ensembl_rel <- "89"

# The genome build used by the Ensembl release
ensembl_genome <- "GRCm38"

# The name of the species used in the filename
ensembl_file_species <- "Mus_musculus"

# The name of the species used in the directory name
ensembl_dir_species <- tolower(ensembl_file_species)

# Setup ------------------------------------------------------------------------

library("Rsubread")

dir.create(outdir, showWarnings = FALSE)

# Construct the URL to the Ensembl FTP site
ensembl_ftp <- paste0("ftp://ftp.ensembl.org/pub/release-", ensembl_rel,
                      "/fasta/", ensembl_dir_species, "/dna/")

# Download chromosome fasta files ----------------------------------------------

for (chr in chroms) {
  cat("Downloading chromesome", chr, "...\n")
  chr_url <- paste0(ensembl_ftp, "/", ensembl_file_species, ".",
                    ensembl_genome, ".dna_sm.chromosome.",
                    chr, ".fa.gz")
  chr_file <- paste0(outdir, "/", ensembl_file_species, ".",
                     ensembl_genome, ".dna_sm.chromosome.",
                     chr, ".fa.gz")
  download.file(url = chr_url, destfile = chr_file)
}

# Decompress and combine fasta files -------------------------------------------

# Name of temporary fasta file
fa_tmp <- paste0("/tmp/tmp.fa")

cmd <- sprintf("zcat %s/*fa.gz > %s", outdir, fa_tmp)
system(cmd)

# Index genome for mappig with Subread -----------------------------------------

buildindex(basename = paste0(outdir, "/", ensembl_genome),
           reference = fa_tmp)

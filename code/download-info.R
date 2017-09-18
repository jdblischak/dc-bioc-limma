#!/usr/bin/env Rscript

# Download SRA information for RNA-seq data from Hoban et al., 2016.
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75816
#
# The output file is GSE75816.txt.
#
# Warning: This downloads a 31 GB file that can be deleted afterwards.
#
# Usage:
#
# Rscript download-info.R outdir
#
# outdir - Directory to save files, created if it doesn't exist.
#
# Steps:
#
# 1. Download SRA database.
#
# 2. Download sample information from GEO.
#
# 3. Obtain URLs to download SRA files from FTP site.
#
# Results:
#
# Writes tab-delimited file GSE75816.txt to outdir.
#
# Also saves the 31 GB SRA sqlite database to outdir.

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("SRAdb"))
suppressPackageStartupMessages(library("stringr"))

# Input command line parameters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript download-info.R outdir")
}
outdir <- args[1]
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# For interactive testing:
# outdir <- "/project2/gilad/jdblischak/dc/data"
# sra_db <- "/project2/gilad/jdblischak/dc/data/SRAmetadb.sqlite"

# 1. Download SRA database and obtain sample information -----------------------

# Download SRAmetadb.sqlite.gz to output directory
sra_db <- file.path(outdir, "SRAmetadb.sqlite")
if (!file.exists(sra_db)) {
  getSRAdbFile(destdir = outdir, method = "wget")
}
stopifnot(file.exists(sra_db))

# 2. Download sample information from GEO --------------------------------------

options('download.file.method.GEOquery' = 'wget')
gse_raw <- getGEO("GSE75816")
gse <- as.data.frame(gse_raw[[1]])
sra_info <- gse %>%
  mutate(id = str_split_fixed(gse$title, " ", 4)[, 1] %>% tolower,
         srx = basename(as.character(supplementary_file_1))) %>%
  select(id, srx)
stopifnot(length(unique(sra_info$id)) == length(sra_info$id),
          length(unique(sra_info$srx)) == length(sra_info$srx))

# 3. Obtain URLs to download SRA files from FTP site ---------------------------

# Translate from SRX to SRR
sra_con <- dbConnect(SQLite(), sra_db)
original_name <- character(length = nrow(sra_info))
new_name <- character(length = nrow(sra_info))
ftp <- character(length = nrow(sra_info))
for (srx_index in 1:nrow(sra_info)) {
  srr <- sraConvert(sra_info$srx[srx_index], "run", sra_con)$run
  # They uploaded one run (i.e. fastq) per sample
  stopifnot(length(srr) == 1)
  sra_ftp_info <- listSRAfile(srr, sra_con)
  ftp[srx_index] <- sra_ftp_info$ftp
  original_name[srx_index] <- paste0(srr, ".sra")
  new_name[srx_index] <- paste0(sra_info$id[srx_index], ".sra")
}

dbDisconnect(sra_con)

# Save SRA information

write.table(data.frame(original_name, new_name, ftp),
            file = file.path(outdir, "GSE75816.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

#!/bin/env Rscript

# IMPORT LIBRARIES
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))

	# DEFINE VARIABLES
reads <- read.csv(file=snakemake@input[["reads"]], sep="\t", header = TRUE)
coverage <- read.csv(file=snakemake@input[["coverage"]], sep="\t", header = TRUE)                                                                                   
# RUN FILTER
source(paste0(snakemake@params[["script_dir"]], "filter_function.R"))
      outputdf <- filter_by_coverage(reads, coverage)
# REMOVE LARGE OBJECTS
rm(reads, coverage)
# SAVE FILTERED CSV
outputdf$Call_ID <- str_replace(outputdf$Call_ID, "_ExpansionHunter", "")
write.table(outputdf, snakemake@output[["csv"]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")



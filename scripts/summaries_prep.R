#!/usr/bin/env Rscript

#IMPORT LIBRARIES
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))

#IMPORT FUNCTIONS
source(paste0(snakemake@params[["script_dir"]], "summaries_functions.R"))

#IMPORT DATA
reads <- read.csv(file=snakemake@input[["csv"]], sep="\t", header = TRUE)
names(reads) <- c("Call_ID", "Sample_ID", "Chr", "Start", "End", "GT", "Ref_Units", "All1", "All2")
outpath <- paste0(dirname(snakemake@output[["csv"]]), "/", snakemake@params[["cohort"]])
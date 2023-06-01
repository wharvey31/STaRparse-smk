#!/usr/bin/env Rscript

#       IMPORT LIBRARIES
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))


    #IMPORT FUNCTIONS
source(paste0(snakemake@params[["script_dir"]], "summaries_functions.R"))

#IMPORT DATA
reads <- read.csv(file=snakemake@input[["reads"]], sep="\t", header = TRUE)
names(reads) <- c("Call_ID", "Sample_ID", "Chr", "Start", "End", "GT", "Ref_Units", "All1", "All2")
outpath <- paste0(snakemake@output[["csv"]]

#RUN FUNCTIONS
locusdf <- by_locus(reads, opt$out, snakemake@wildcards[["sample"]], opt$build)
ratio_of_stability(locusdf, outpath)
by_chr(reads, locusdf, outpath)
by_sample(reads, locusdf, outpath)

#!/usr/bin/env Rscript

#       IMPORT LIBRARIES
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))

#       DEFINE ARGUMENTS FOR SCRIPT
option_list <- list(
  make_option(c("-d", "--dir"), type="character", action="store", default=NA, help="Path of global directory", metavar="character"),
  make_option(c("-r", "--reads"), type="character", action="store", default=NA, help="Read file path", metavar="character"),
  make_option(c("-o", "--out"), type="character", action="store", default=NA, help="output file path and name", metavar="character"),
  make_option(c("-s", "--savename"), type="character", action="store", help="Indicate job name", metavar="character"),
  make_option(c("-b", "--build"), type="character", action="store", help="Specify reference genome build (19, 37, or 38)", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

#       ARGUMENT CHECK

  
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

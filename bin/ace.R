#!/usr/bin/env Rscript

library(ACE)
args <- commandArgs(trailingOnly = TRUE)
output_folder <- args[1]
genome <- if (length(args) >= 2) args[2] else "hg38"
ploidy_arg <- if (length(args) >= 3) args[3] else "2,3,4"
ploidy <- as.numeric(strsplit(ploidy_arg, ",", fixed = TRUE)[[1]])

if (any(is.na(ploidy))) {
  stop(sprintf("Invalid ACE ploidy value '%s'. Use a comma-separated numeric list, e.g. 2 or 2,3,4.", ploidy_arg))
}


runACE(
  outputdir = output_folder, filetype = "bam",
  binsizes = 1000, genome = genome, ploidies = ploidy
)

# runACE(inputdir = "./", outputdir, filetype = 'rds', genome = 'hg19',
#  binsizes, ploidies = 2, imagetype = 'pdf', method = 'RMSE', penalty = 0,
#  cap = 12, bottom = 0, trncname = FALSE, printsummaries = TRUE,
#  autopick = FALSE)

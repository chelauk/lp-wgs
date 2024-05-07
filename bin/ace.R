#!/usr/bin/env Rscript

library(ACE)
args <- commandArgs(trailingOnly = TRUE)
output_folder <- args[1]
ploidy <- c(2, 3, 4)


runACE(
    outputdir = output_folder, filetype = "bam",
    binsizes = 1000, genome = "hg38", ploidies = ploidy
)

# runACE(inputdir = "./", outputdir, filetype = 'rds', genome = 'hg19',
#  binsizes, ploidies = 2, imagetype = 'pdf', method = 'RMSE', penalty = 0,
#  cap = 12, bottom = 0, trncname = FALSE, printsummaries = TRUE,
#  autopick = FALSE)

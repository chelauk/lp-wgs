#!/usr/bin/env Rscript

library(ACE)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
output_folder <- args[1]
genome <- if (length(args) >= 2) args[2] else "hg38"
ploidy_arg <- if (length(args) >= 3) args[3] else "2,3,4"
binsize <- if (length(args) >= 4) as.numeric(args[4]) else 1000
binsize_label <- format(binsize, trim = TRUE, scientific = FALSE)
inputdir <- if (length(args) >= 5) args[5] else "."
ploidy <- as.numeric(strsplit(ploidy_arg, ",", fixed = TRUE)[[1]])

if (any(is.na(ploidy))) {
  stop(sprintf("Invalid ACE ploidy value '%s'. Use a comma-separated numeric list, e.g. 2 or 2,3,4.", ploidy_arg))
}

if (is.na(binsize)) {
  stop(sprintf("Invalid ACE bin size '%s'. Use a numeric QDNAseq bin size, e.g. 1000.", args[4]))
}

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

rds_pattern <- paste0(".*", binsize_label, "kbp\\.rds$")
rds_files <- list.files(inputdir, pattern = rds_pattern, full.names = TRUE)

if (length(rds_files) != 1) {
  stop(sprintf("Expected exactly one QDNAseq RDS matching '%s' in '%s', found %s.",
               rds_pattern, inputdir, length(rds_files)))
}

qdnaseq_object <- readRDS(rds_files[1])
sample_prefix <- sub(paste0("_", binsize_label, "kbp\\.rds$"), "",
                     basename(rds_files[1]))

# Use ACE to estimate minima from the prepared QDNAseq object. Full runACE()
# below writes the standard ACE output directory for downstream consumers.
template <- objectsampletotemplate(qdnaseq_object, index = 1)
sqmodel <- squaremodel(template, penalty = 0.5, penploidy = 0.5)
pdf(file.path(output_folder, paste0(sample_prefix, "_sky_on_fire.pdf")))
sqmodel$matrixplot + ggtitle(paste0(sample_prefix, "_Sky on fire"))
dev.off()

write.table(file = file.path(output_folder,
                             paste0(sample_prefix, "_sqmodel_minmadf.txt")),
            quote = FALSE, x = sqmodel$minimadf)

runACE(
  inputdir = inputdir, outputdir = output_folder, filetype = "rds",
  binsizes = binsize, genome = genome, ploidies = ploidy
)

# runACE(inputdir = "./", outputdir, filetype = 'rds', genome = 'hg19',
#  binsizes, ploidies = 2, imagetype = 'pdf', method = 'RMSE', penalty = 0,
#  cap = 12, bottom = 0, trncname = FALSE, printsummaries = TRUE,
#  autopick = FALSE)

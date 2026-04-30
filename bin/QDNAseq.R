#!/usr/bin/env Rscript

# Running QDNAseq for low pass analysis
#
# Coded by George Cresswell; 2018-05-11
# Modified by Chela James;   2024-08-03

args <- commandArgs(trailingOnly = TRUE)

if (!require(QDNAseq)) stop("Package 'QDNAseq' missing\n.")
if (!require(CGHcall)) stop("Package 'CGHcall' missing\n.")
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")
if (!require(ACE)) stop("Package 'ACE' missing\n.")

patient <- args[1]
sample  <- args[2]
binsize <- as.numeric(args[3])
bam     <- args[4]
genome  <- if (length(args) >= 5) args[5] else "hg38"
qdnaseq_package <- if (length(args) >= 6) args[6] else paste0("QDNAseq.", genome)

if (!require(qdnaseq_package, character.only = TRUE)) {
  stop(sprintf("Package '%s' missing\n.", qdnaseq_package))
}

autosomes <- switch(
  genome,
  hg38 = as.character(1:22),
  hg19 = as.character(1:22),
  mm10 = as.character(1:19),
  stop(sprintf("Unsupported QDNAseq genome '%s'", genome))
)

normalize_chr <- function(x) {
  gsub("^chr", "", as.character(x), ignore.case = TRUE)
}

bins <- getBinAnnotations(binSize = binsize, genome = genome)

# These steps are taken directly from the QDNAseq tutorial
read_counts <- binReadCounts(bins, bamfiles = bam)
read_counts_filtered <- applyFilters(read_counts, residual = TRUE,
                                     blacklist = TRUE, mappability = 65,
                                     bases = 95)
read_counts_filtered <- estimateCorrection(read_counts_filtered)
read_counts_filtered_xy <- applyFilters(read_counts_filtered, residual = TRUE,
                                        blacklist = TRUE, mappability = 65,
                                        bases = 95, chromosomes = NA)
copy_numbers <- correctBins(read_counts_filtered_xy)
copy_numbers_normalized <- normalizeBins(copy_numbers)

# copy_numbers_smooth <- smoothOutlierBins(copy_numbers_normalized), # nolint
# removed as it may smooth real amplifications

# Has a random element, so set seed
set.seed(1)
# Segment the data
copy_numbers_segmented <- segmentBins(copy_numbers_normalized,
                                      transformFun = "sqrt")

# SegmentBins implements the DNAcopy::segment function which uses the
# circular binary segmentation (CBS) algorithm. This algorithm is not
# deterministic, so you may get different results each time you run it.
# Setting the seed ensures that you get the same results each time you
# run the code.

# in short it splits the chromosome into segments of similar copy number,
# and then estimates the mean copy number for each segment.

# 0.01  0.03 -0.02  0.00 | 0.55  0.61  0.58  0.60

# Low-count bins and high-count bins do not have exactly the same variance.
# A square-root transform is a classic variance-stabilising transform for count-like data.

copy_numbers_segmented <- normalizeSegmentedBins(copy_numbers_segmented)

# use the sqmodel function to estimate minima
# You'll find that ACE can make great fits at interesting ploidies,
# but they may not always make sense from a biological perspective.
# On top of the penalty for low cellularities, you can consider to
# use a penalty for ploidies (penploidy) that diverge a lot from two.
# To use data from QDNAseq-objects, ACE parses it into data frames
# referred to as "templates". Because we will look at sample2 several
# times, we can just create a variable with this data frame.

template <- objectsampletotemplate(copy_numbers_segmented, index = 1)
sqmodel <- squaremodel(template, penalty = 0.5, penploidy = 0.5)
pdf(paste0(patient, "_", sample, "_sky_on_fire.pdf"))
sqmodel$matrixplot + ggtitle(paste0(patient, "_", sample, "_Sky on fire"))
dev.off()

write.table(file = paste0(patient, "_", sample, "_sqmodel_minmadf.txt"),
            quote = FALSE, x = sqmodel$minimadf)

copy_numbers_called <- callBins(copy_numbers_segmented, method = c("cutoff"),
                                cutoffs = log2(c(deletion = 2 - 1.2,
                                                 loss = 2 - 0.2, gain = 2 + 0.2,
                                                 amplification = 2 + 1.2) / 2))

# Out the plot of called segments
pdf(paste0(patient, "_", sample, "_called_segments.pdf"))
plot(copy_numbers_called, ylim = c(-2, 2))
abline(v = log2(c(deletion = 2 - 1.2, loss = 2 - 0.2, gain = 2 + 0.2,
                  amplification = 2 + 1.2) / 2))
dev.off()

# outputBINS
print(paste("bins_file:", paste0(patient, "_", sample, "_bins.txt")))
exportBins(copy_numbers_called,
           file = paste0(patient, "_", sample, "_bins.txt"))

# Convert to a CGHcall style object
cgh <- makeCgh(copy_numbers_called)

# Add two to calls to make them like they are diploid, for comparison
calls_adj <- calls(cgh) + 2

# Write out the calls for future fun
write.table(calls_adj, file = paste0(patient, "_", sample, "_cna_calls.txt"),
            quote = FALSE, sep = "\t")

# Write out the segments with mean logr for extra fun
write.table(segmented(cgh),
            file = paste0(patient, "_", sample, "_cna_segments.txt"),
            quote = FALSE, sep = "\t")

# Read in the data
lrrs <- read.table(paste0(patient, "_", sample, "_bins.txt"), header = TRUE)
call <- read.table(paste0(patient, "_", sample, "_cna_calls.txt"),
                   header = TRUE)

# Play with levels
call <- as.factor(call[, 1])
call <- factor(call, levels = c("-2", "-1", "0", "1", "2"))

segs <- read.table(paste0(patient, "_", sample, "_cna_segments.txt"),
                   header = TRUE)
med_lrrs <- median(lrrs[normalize_chr(lrrs$chromosome) %in% autosomes, 5])

# Adjust the lrrs _and_ the segment lrrs by this value
lrrs[, 5] <- lrrs[, 5] - med_lrrs
segs      <- segs - med_lrrs

# Return them back outside, weird way of doing it, but it's because the
# original export is embedded in a QDNAseq function
write.table(lrrs, file = paste0(patient, "_", sample, "_bins.txt"),
            quote = FALSE, sep = "\t")
write.table(segs, file = paste0(patient, "_", sample, "_cna_segments.txt"),
            quote = FALSE, sep = "\t")

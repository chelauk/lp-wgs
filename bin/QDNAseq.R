#!/usr/bin/env Rscript

#!Rscript

# Running QDNAseq for low pass analysis
# 
# Coded by George Cresswell; 2018-05-11
# Modified by Chela James;   2024-08-03

args = commandArgs(trailingOnly=TRUE)

if (!require(QDNAseq)) stop("Package 'QDNAseq' missing\n.")
if (!require(CGHcall)) stop("Package 'CGHcall' missing\n.")
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")
if (!require(QDNAseq.hg38)) stop("Package 'QDNAseq.hg38' missing\n.")
if (!require(ACE)) stop("Package 'ACE' missing\n.")

patient = args[1]
sample  = args[2]
binsize = as.numeric(args[3])
bam     = args[4]

bins <- getBinAnnotations(binSize = binsize, genome = "hg38")

# These steps are taken directly from the QDNAseq tutorial
readCounts = binReadCounts(bins, bamfiles=bam)
readCountsFiltered = applyFilters(readCounts, residual=TRUE, blacklist=TRUE, mappability = 65, bases = 95)
readCountsFiltered = estimateCorrection(readCountsFiltered)
readCountsFiltered_XY = applyFilters(readCountsFiltered, residual=TRUE, blacklist=TRUE, mappability = 65, bases = 95, chromosomes=NA)
copyNumbers = correctBins(readCountsFiltered_XY)
copyNumbersNormalized = normalizeBins(copyNumbers)
# copyNumbersSmooth = smoothOutlierBins(copyNumbersNormalized), removed as it may smooth real amplifications

# Has a random element, so set seed
set.seed(1) 
# Segment the data
copyNumbersSegmented = segmentBins(copyNumbersNormalized, transformFun="sqrt")
copyNumbersSegmented = normalizeSegmentedBins(copyNumbersSegmented)

# use the sqmodel function to estimate minima
# You'll find that ACE can make great fits at interesting ploidies, but they may not always make sense from
# a biological perspective. 
# On top of the penalty for low cellularities, you can consider to use a penalty for ploidies (penploidy) 
# that diverge a lot from two.
# To use data from QDNAseq-objects, ACE parses it into data frames referred to
# as "templates". Because we will look at sample2 several times, we can just
# create a variable with this data frame.

template <- objectsampletotemplate(copyNumbersSegmented, index = 1)
sqmodel <- squaremodel(template, penalty = 0.5, penploidy = 0.5,)
pdf(paste0(patient,"_",sample,"_sky_on_fire.pdf"))
sqmodel$matrixplot + ggtitle(paste0(patient,"_",sample,"_Sky on fire"))
dev.off()

write.table(file = paste0(patient,"_",sample,"_sqmodel_minmadf.txt") , quote = FALSE, x = sqmodel$minimadf)

copyNumbersCalled = callBins(copyNumbersSegmented,
                             method=c("cutoff"),
                             cutoffs=log2(c(deletion = 2 - 1.2, loss = 2 - 0.2, gain = 2 + 0.2, amplification = 2 + 1.2)/2))

# Out the plot of called segments
pdf(paste0("EX24_",sample,"_called_segments.pdf"))
  plot(copyNumbersCalled, ylim = c(-2,2))
  abline(v=log2(c(deletion = 2 - 1.2, loss = 2 - 0.2, gain = 2 + 0.2, amplification = 2 + 1.2)/2))
dev.off()

# outputBINS
print(paste("bins_file:",paste0("EX24_",sample,"_bins.txt")))
exportBins(copyNumbersCalled, file = paste0("EX24_",sample,"_bins.txt"))


# Convert to a CGHcall style object
cgh = makeCgh(copyNumbersCalled)

# Add two to calls to make them like they are diploid, for comparison
calls_adj = calls(cgh)+2

# Write out the calls for future fun
write.table(calls_adj, file = paste0("EX24_",sample,"_cna_calls.txt"), quote = FALSE, sep = "\t")

# Write out the segments with mean logr for extra fun
write.table(segmented(cgh), file = paste0("EX24_",sample,"_cna_segments.txt"), quote = FALSE, sep = "\t")

# Read in the data
lrrs = read.table(paste0("EX24_",sample,"_bins.txt"), header = TRUE)
call = read.table(paste0("EX24_",sample,"_cna_calls.txt"), header = TRUE)

# Play with levels
call = as.factor(call[,1])
call = factor(call, levels = c("-2", "-1", "0", "1", "2"))

segs = read.table(paste0("EX24_",sample,"_cna_segments.txt"), header = TRUE)
med_lrrs = median(lrrs[lrrs$chromosome %in% 1:22,5])

# Adjust the lrrs _and_ the segment lrrs by this value
lrrs[,5] = lrrs[,5] - med_lrrs
segs     = segs - med_lrrs

# Return them back outside, weird way of doing it, but it's because the original export is embedded in a QDNAseq function
write.table(lrrs, file = paste0("EX24_",sample,"_bins.txt"), quote = FALSE, sep = "\t")
write.table(segs, file = paste0("EX24_",sample,"_cna_segments.txt"), quote = FALSE, sep = "\t")
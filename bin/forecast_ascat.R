library(copynumber)
library(ggplot2)
library(cowplot)
source("scripts/00_general_functions.R")
source("scripts/runASCATlp.R")

args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]
samples <- unlist(strsplit(args[2], " "))
ids     <- unlist(strsplit(args[3], " "))
ploidy  <- as.numeric(args[4])
purity  <- as.numeric(args[5])

source(paste0(bin_dir, "/00_general_functions.R"))
source(paste0(bin_dir, "/runASCATlp.R"))

# These are the arms of hg38
arms <- read.table(bin_dir, "/chrArmBoundaries_hg38.txt", header = TRUE)

# log2R data output by QDNAseq.R
# https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf
# bins <- getBinAnnotations(binSize = binsize, genome = "hg38")
# readCounts = binReadCounts(bins, bamfiles=bam)
# readCountsFiltered = applyFilters(readCounts, residual=TRUE, blacklist=TRUE, mappability = 65, bases = 95)
# readCountsFiltered = estimateCorrection(readCountsFiltered)
# readCountsFiltered_XY = applyFilters(readCountsFiltered, residual=TRUE, blacklist=TRUE, mappability = 65, bases = 95, chromosomes=NA)
# copyNumbers = correctBins(readCountsFiltered_XY)
# copyNumbersNormalized = normalizeBins(copyNumbers)
# copyNumbersSegmented = segmentBins(copyNumbersNormalized, transformFun="sqrt")
# copyNumbersSegmented = normalizeSegmentedBins(copyNumbersSegmented)
# copyNumbersCalled = callBins(copyNumbersSegmented,
#                              method=c("cutoff"),
#                              cutoffs=log2(c(deletion = 2 - 1.2, loss = 2 - 0.2, gain = 2 + 0.2, amplification = 2 + 1.2)/2))

## outputBINS
# print(paste("bins_file:",paste0(patient,"_",sample,"_bins.txt")))
# exportBins(copyNumbersCalled, file = paste0(patient,"_",sample,"_bins.txt"))
# lrrs = read.table(paste0(patient, "_",sample,"_bins.txt"), header = TRUE)
# med_lrrs = median(lrrs[lrrs$chromosome %in% 1:22,5])
# lrrs[,5] = lrrs[,5] - med_lrrs
# write.table(lrrs, file = paste0(patient,"_",sample,"_bins.txt"), quote = FALSE, sep = "\t")


patient_lrr <- read.table(paste0(id, "_bins.txt"),
                          header = TRUE, stringsAsFactors = FALSE)

# Get the index of bins
chr_pos <- patient_lrr[, 2:4]
bin_row <- patient_lrr[, 1]

# chr_pos
chr_pos$chromosome <- factor(chr_pos$chromosome, levels = c(1:22, "X", "Y"))

# Take only the lrr data
p_mat <- cbind(patient_lrr[, 5])
colnames(p_mat) <- id

# function pcf = Single-sample copy number segmentation
# gamma
# penalty for each discontinuity in the curve, default is 40.
sample_seg <- pcf(data.frame(chr = chr_pos$chromosome,
                             pos = chr_pos$start,
                             p_mat),
                  arms = getPQ(chr_pos, arms),
                  gamma = 10, fast = FALSE)

colnames(sample_seg)[7] <- sample_seg$sampleID[1]
sample_seg$sampleID <- NULL

# now we expand according to the number of probes
expanded_segs <- rep(sample_seg[, 6], times = sample_seg$n.probes)
names(expanded_segs) <- rep(id, length(expanded_segs))

# runAscat
mid_pld <- 3.1
expand  <- 1.6
mp      <- 1

autosome_index <- chr_pos$chromosome %in% 1:22
bins_auto <- p_mat[autosome_index]
segs_auto <- expanded_segs[autosome_index]
sn <- id
ps <- FALSE
pr <- 1000
pp <- 1000

res <- runASCATlp(lrrs = segs_auto, fix_ploidy = mid_pld, pad_ploidy = expand,
                  interval = 0.01, min_purity = 0.01, max_lrr = Inf,
                  no_fit_psit = 2, preset = ps, preset_purity = pr,
                  preset_ploidy = pp, max_purity = mp)


sex_segs <- expanded_segs[!autosome_index]
n <- callXchromsome(sex_lrrs = sex_segs, psi = res$Psi, psit = res$PsiT,
                    purity = res$Purity)

res$CN <- c(res$CN, ifelse(round(n) < 0, 0, round(n)))
res$contCN <- c(res$contCN, n)
res$segs <- expanded_segs
res$bins <- patient_lrr[, 5]
res$sample_name <- id

cna_data    <- res

sample_name <- cna_data$samplename

cn_output   <- data.frame(cna_data$CN)
colnames(cn_output) <- sample_name

# Make a plot dataframe
plt_df <- data.frame(genome.bin = seq_along(cna_data$bins),
                     Chromosome = chr_pos$chromosome,
                     Log2ratio = cna_data$bins,
                     mean_segment = cna_data$segs,
                     Call = as.factor(cna_data$CN))

# Max CN
max_cn <- 2
min_cn <- -2

# lines across
lines_across <- data.frame(x1 = 1,
                           y1 = min_cn:max_cn,
                           x2 = length(cna_data$bins),
                           y2 = min_cn:max_cn)

# lines going up for chromosomes
lines_vertical <- data.frame(x1 = c(1, cumsum(table(plt_df$Chromosome))),
                             y1 = min_cn,
                             x2 = c(1, cumsum(table(plt_df$Chromosome))),
                             y2 = max_cn)

# Remove any chromosome labels due to congestion?
chr_out <- c(19, 21)

chrs_lab <- c(1:22, "X", "Y")
chrs_lab[chr_out] <- ""

# Make the plot
p <- ggplot(plt_df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
  geom_hline(yintercept = c(-2, -1, 0, 1, 2), lty = c("solid"), lwd = 0.2) +
  geom_point() +
  scale_colour_manual(values = cols) +
  scale_x_continuous(name = "Chromosomes", labels = chrs_lab,
                     breaks = as.vector(c(1, cumsum(table(plt_df$Chromosome))[-24]) + # nolint: line_length_linter.
                                          (table(plt_df$Chromosome) / 2))) +
  geom_vline(xintercept = c(1, cumsum(table(plt_df$Chromosome))),
             lty = "dotted") +
  ggtitle(paste0("Low pass calls - ", sample_name, ", purity=",
                 cna_data$Purity, ", psit = ", cna_data$PsiT)) +
  scale_y_continuous(limits = c(-2, 2), oob = scales::squish) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18)) +
  geom_point(aes(y = mean_segment), color = "#000000")

autosome_index <- chr_pos$chromosome %in% 1:22

ploidy_expected <- rep(1, times = length(cna_data$CN))
ploidy_expected[!autosome_index] <- 0.5
ploidy_expected <- median(cna_data$CN) * ploidy_expected
ploidy_expected <- round(ploidy_expected)

# Get output metrics
metrics <- data.frame(Sample = sample_name, Purity = cna_data$Purity,
                      PsiT = cna_data$PsiT,
                      Ploidy = signif(mean(cna_data$CN[autosome_index]),
                                      digits = 4),
                      PGA = signif(length(which(
                                                cna_data$CN != ploidy_expected))
                      / length(cna_data$CN),
                      digits = 4))


cn_out <- data.frame(res$CN)

colnames(cn_out) <- names(ascat)[i]
rownames(cn_out) <- patient_lrr$feature
write.table(cn_out, file = paste0(id, "_cna_ploidy_search_calls.txt"),
            quote = FALSE)
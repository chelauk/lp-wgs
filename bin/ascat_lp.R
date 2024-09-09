#!/usr/bin/env Rscript

library(copynumber)
library(ggplot2)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]
samples <- unlist(strsplit(args[2]," "))
ids     <- unlist(strsplit(args[3]," "))
ploidy  <- as.numeric(args[4])
purity  <- as.numeric(args[5])
bin_dir <- args[6]

source(paste0(bin_dir,"/00_general_functions.R"))
source(paste0(bin_dir,"/runASCATlp.R"))

# create output dir
if (!dir.exists("./plots")) {
  # Create the directory if it does not exist
  dir.create("./plots")
  message("Directory created: ", "./plots")
} else {
  message("Directory already exists: ", "./plots")
}

# GRCh38 arms
arms <- read.table(paste0(bin_dir,"/chrArmBoundaries_hg38.txt"))

files <- list.files('.', pattern = "_bins.txt")

# Read the files in and store as list
per_patient_lrrs = lapply(files, function(fi) {

  lapply(fi, function(f) {
    # Extract the unique part of the filename to use as the last column name
    file_id <- gsub("_bins\\.txt$", "", basename(f))
    # Define the column names
    col_names <- c("feature", "chromosome", "start", "end", file_id)
    read.table(f, stringsAsFactors = F, col.names = col_names)
  })
})

# Get the index of bins
chr_pos = per_patient_lrrs[[1]][[1]][,2:4]
bin_row = per_patient_lrrs[[1]][[1]][,1]

# chr_pos
chr_pos$chromosome = factor(chr_pos$chromosome, levels = c(1:22,"X","Y"))

# Take only the lrr data
per_patient_lrrs = lapply(per_patient_lrrs, function(p) {

  lapply(p, function(l) {

    cn = colnames(l)[5]

    l = cbind(l[,5])

    colnames(l) = cn

    return(l)

  })

})

# Make it a matrix
per_patient_lrr_mat = lapply(per_patient_lrrs, function(p) do.call(cbind, p))

# Perform multi-sample segmentation calling
per_patient_ms_segs = lapply(per_patient_lrr_mat, function(p_mat) {

  if(ncol(p_mat) > 1) {

    # Perform the multiregion segmentation
    per_patient_multisample_seg = multipcf(data.frame(chr = chr_pos$chromosome,
                                                      pos = chr_pos$start,
                                                      p_mat),
                                           arms = getPQ(chr_pos, arms),
                                           gamma = 10, fast = FALSE)



  } else {

    per_patient_multisample_seg = pcf(data.frame(chr = chr_pos$chromosome,
                                                 pos = chr_pos$start,
                                                 p_mat),
                                      arms = getPQ(chr_pos, arms),
                                      gamma = 10, fast = FALSE)

    colnames(per_patient_multisample_seg)[7] = per_patient_multisample_seg$sampleID[1]
    per_patient_multisample_seg$sampleID = NULL

  }

  # Expand it out per sample
  per_patient_ms_segs = lapply(1:(ncol(per_patient_multisample_seg) - 5), function(s) {

    rep(per_patient_multisample_seg[,s+5], times = per_patient_multisample_seg[,5])

  })

  # Name the elements in the list
  names(per_patient_ms_segs) = colnames(per_patient_multisample_seg)[6:ncol(per_patient_multisample_seg)]

  return(per_patient_ms_segs)

})

# Run lpASCAT
per_patient_ms_ascat = lapply(1:length(per_patient_ms_segs), function(i) {

  p = per_patient_ms_segs[[i]]

  lapply(1:length(p), function(s) {

    # Write out the segments for use in the manual ploidy choosing script
    samples = names(p)

    mid_pld = 3.1
    expand  = 1.6
    mp      = 1

    if(any(grepl("XXXXX", samples))) {

      mid_pld = 4.35
      expand  = 0.35

    }

    autosome_index = chr_pos$chromosome %in% 1:22

    bins_auto = per_patient_lrr_mat[[i]][autosome_index,s]
    segs_auto = p[[s]][autosome_index]

    sn = samples[s]

    ps = F
    pr = 1000
    pp = 1000

    if(sn == "XXXXX") {

      ps = T
      pr = 0.07
      pp = 4.41

    }

    res = runASCATlp(lrrs = segs_auto, fix_ploidy = mid_pld, pad_ploidy = expand,
                     interval = 0.01, min_purity=0.01, max_lrr = Inf, no_fit_psit = 2, preset = ps,
                     preset_purity = pr, preset_ploidy = pp, max_purity = mp)

    sex_segs = p[[s]][!autosome_index]

    n = callXchromsome(sex_lrrs = sex_segs, psi = res$Psi, psit = res$PsiT, purity = res$Purity)

    res$CN = c(res$CN, ifelse(round(n) < 0, 0, round(n)))
    res$contCN = c(res$contCN, n)

    res$segs       = p[[s]]
    res$bins       = per_patient_lrr_mat[[i]][,s]
    res$samplename = sn

    return(res)

  })

})

# Make plots per sample
lapply(per_patient_ms_ascat, function(i) {

  lapply(i, function(j) {

    cna_data = j

    sample_name = cna_data$samplename
    #print(sample_name)

    cn_output   = data.frame(cna_data$CN)
    colnames(cn_output) = sample_name

    # Make a plot dataframe
    plt.df = data.frame(genome.bin = 1:length(cna_data$bins),
                        Chromosome = chr_pos$chromosome,
                        Log2ratio = cna_data$bins,
                        mean_segment = cna_data$segs,
                        Call = as.factor(cna_data$CN))
    # Max CN
    maxCN = 2
    minCN = -2

    # lines across
    lines_across = data.frame(x1 = 1,
                              y1 = minCN:maxCN,
                              x2 = length(cna_data$bins),
                              y2 = minCN:maxCN)

    # lines going up for chromosomes
    lines_vertical = data.frame(x1 = c(1, cumsum(table(plt.df$Chromosome))),
                                y1 = minCN,
                                x2 = c(1, cumsum(table(plt.df$Chromosome))),
                                y2 = maxCN)

    # Remove any chromosome labels due to congestion?
    chr_out = c(19,21)

    chrs_lab = c(1:22,"X","Y")
    chrs_lab[chr_out] = ""
    print(paste("title:",paste0("Low pass calls - ",sample_name,", purity=",
                                cna_data$Purity,", psit = ",cna_data$PsiT)))
    # Make the plot
    p = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_hline(yintercept = c(-2,-1,0,1,2), lty = c("solid"), lwd = 0.2) +
      geom_point() +
      scale_colour_manual(values = cols) +
      scale_x_continuous(name = "Chromosomes", labels = c(1:22,"X","Y"),
                         breaks = as.vector(c(1, cumsum(table(plt.df$Chromosome))[-24]) +
                                              (table(plt.df$Chromosome) / 2))) +
      geom_vline(xintercept = c(1, cumsum(table(plt.df$Chromosome))), lty = "dotted") +
      ggtitle(paste0("Low pass calls - ",sample_name,", purity=",
                     cna_data$Purity,", psit = ",cna_data$PsiT)) +
      scale_y_continuous(limits=c(-2,2), oob=scales::squish) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18)) +
      geom_point(aes(y = mean_segment), color="#000000")
    ggsave(filename = paste0("./plots/",sample_name,"_multiregion_seg_cna_profile.png"),
           plot = p, width = 15, height = 5)

    write.table(cn_output, file = paste0("./plots/",sample_name,"_multiregion_seg_cna_calls.txt"),
                row.names = T, col.names = T, quote = F)

    # Calculate ploidy vector
    autosome_index = chr_pos$chromosome %in% 1:22

    ploidy_expected = rep(1, times = length(cna_data$CN))
    ploidy_expected[!autosome_index] = 0.5
    ploidy_expected = median(cna_data$CN)*ploidy_expected
    ploidy_expected = round(ploidy_expected)

    # Get output metrics
    metrics = data.frame(Sample = sample_name, Purity = cna_data$Purity, PsiT = cna_data$PsiT,
                         Ploidy = signif(mean(cna_data$CN[autosome_index]), digits = 4),
                         PGA = signif(length(which(cna_data$CN!=ploidy_expected)) / length(cna_data$CN), digits = 4))

    if(median(cna_data$CN) != round(metrics$Ploidy)) {warnings("Median CN and mean ploidy mismatch!")}

    write.table(metrics, file = paste0("./plots/",sample_name,"_multiregion_seg_metrics.txt"),
                row.names = F, col.names = T, quote = F, sep = "\t")

    # Make plot
    p2 = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_across, color="#00000080", lwd = 0.2) +
      scale_colour_manual(values = cols) +
      scale_x_continuous(name = NULL,
                         labels = chrs_lab,
                         breaks = as.vector(c(1, cumsum(table(plt.df$Chromosome))[-24]) +
                                              (table(plt.df$Chromosome) / 2)),
                         expand = c(0.01,0)) +
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_vertical, color="#0000001A", lwd = 0.2) +
      scale_y_continuous(limits=c(-2,2), oob=scales::squish, expand = c(0,0)) +
      geom_point(aes(size = 2)) +
      geom_point(aes(y = mean_segment), color="#000000") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x = element_text(size = 18),
            #axis.text.x = element_blank(),
            axis.text.x = element_text(size = 40, angle = 45, vjust = 1, hjust=1),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 60),
            legend.position = "none",
            plot.margin = margin(0.4, 0.1, 0.4, 0.1, "in"))
    pdf(paste0("./plots/",sample_name,"_multiregion_seg_cna_profile_separate.pdf"),
        height = 5, width = 24, useDingbats=FALSE)
    dev.off()

  })

})


#!/usr/bin/env Rscript


if (!require(QDNAseq)) stop("Package 'QDNAseq' missing\n.")
if (!require(CGHcall)) stop("Package 'CGHcall' missing\n.")
if (!require(ACE)) stop("Package 'ACE' missing\n.")
if (!require(bcp)) stop("Package 'bcp' missing\n.")
if (!require(tidyverse)) stop("Package 'tidyverse' missing\n")
if (!require(magrittr)) stop("Package 'magrittr' missing\n.")
if (!require(pracma)) stop("Package 'pracma' missing\n")

args <- commandArgs(trailingOnly = TRUE)

patient <- args[1]
sample  <- args[2]
binsize <- as.numeric(args[3])
bin_dir <- args[4]
bam     <- args[5]
genome  <- if (length(args) >= 6) args[6] else "hg38"
qdnaseq_package <- if (length(args) >= 7) args[7] else paste0("QDNAseq.", genome)

if (!require(qdnaseq_package, character.only = TRUE)) {
  stop(sprintf("Package '%s' missing\n.", qdnaseq_package))
}

source(paste0(bin_dir,"/segmentation.R"))
source(paste0(bin_dir,"/helper_functions.R"))
source(paste0(bin_dir,"/00_general_functions.R"))
source(paste0(bin_dir,"/runASCATlp.R"))

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


# create a bins dataframe

bins_df <- data.frame(
  chromosome = normalize_chr(chromosomes(copy_numbers_normalized)),
  start      = bpstart(copy_numbers_normalized),
  end        = bpend(copy_numbers_normalized),
  use        = copy_numbers_normalized@featureData@data$use
)

dir.create(paste0(patient, "_", sample, "_bcp"), showWarnings = FALSE)

for (chr in autosomes) {

  chr_data <- get_chr_signal(copy_numbers_normalized, chr)

  denoised <- bcp(
    chr_data$y,
    p0 = 0.01,
    burnin = 500,
    mcmc = 2000,
    return.mcmc = TRUE
  )

  pdf(
    file = paste0(
      patient, "_", sample,"_bcp/",
      patient, "_", sample, "_bcp_chr", chr, ".pdf"
    ),
    width = 10,
    height = 6
  )

  plot(
    denoised,
    main = paste0(
      "Posterior Means and Probabilities of a Change — chr", chr
    )
  )

  dev.off()
}

bayes_df <- lapply(autosomes, function(chr) {
  run_bayes_chr(copy_numbers_normalized, chr)
}) |> bind_rows()

chr_sizes <- bayes_df |>
  group_by(chromosome) |>
  summarise(chr_len = max(end), .groups = "drop") |>
  arrange(as.numeric(chromosome)) |>
  mutate(offset = dplyr::lag(cumsum(as.numeric(chr_len)), default = 0))

bayes_df <- bayes_df |>
  left_join(chr_sizes, by = "chromosome") |>
  mutate(
    genome_pos = (start + end) / 2 + offset
  )

axis_df <- bayes_df |>
  group_by(chromosome) |>
  summarise(center = mean(genome_pos), .groups = "drop")

p <- ggplot(bayes_df, aes(x = genome_pos)) +
  geom_point(aes(y = raw), color = "grey80", size = 0.3) +
  geom_line(aes(y = profile), color = "firebrick", linewidth = 0.4) +
  geom_vline(
    data = chr_sizes[-1, ],
    aes(xintercept = offset),
    linetype = "dashed",
    color = "grey70"
  ) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chromosome
  )

# convert bayes_df for ascat lp

bayes_segments <- bayes_df |>
  arrange(as.numeric(chromosome), start) |>
  group_by(chromosome) |>
  mutate(
    seg_id = cumsum(breakpoint) + 1
  ) |>
  group_by(chromosome, seg_id) |>
  summarise(
    start = min(start),
    end = max(end),
    n_bins = n(),
    mean_cn = mean(raw, na.rm = TRUE),
    bayes_cn = mean(profile, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(bayes_segments,
          file = paste0(
                        patient, "_", sample,"_bcp/",
                        patient,"_",sample,"_bcp_segments.csv"
                        ),
          row.names = FALSE)

bayes_segments <- bayes_segments |>
  mutate(lrr = log2(mean_cn / 2))

segs_auto <- rep(bayes_segments$lrr, times = bayes_segments$n_bins)

res <- runASCATlp(
  lrrs = segs_auto,
  fix_ploidy = 2,
  pad_ploidy = 1.6,
  interval = 0.01,
  min_purity = 0.01,
  max_lrr = Inf,
  no_fit_psit = 2,
  preset = FALSE,
  preset_purity = 1000,
  preset_ploidy = 1000,
  max_purity = 1
)

solutions <- res$near_best

if (is.null(solutions) || nrow(solutions) == 0) {
  solutions <- data.frame(
    psit = res$PsiT,
    purity = res$Purity,
    dist = NA
  )
}

for (i in seq_len(nrow(solutions))) {

  solution_name <- ifelse(i == 1, "selected", paste0("alt_", i - 1))

  pdf(
    paste0(
      patient, "_", sample, "_bcp/",
      patient, "_", sample, "_bcp_wgs_profile_", solution_name, ".pdf"
    ),
    width = 10,
    height = 6
  )

  print(
    p +
      labs(
        x = "Chromosome",
        y = "Copy number",
        title = paste0(
          "BayesCNA genome-wide profile; ASCATlp ",
          solution_name,
          ": purity = ", solutions$purity[i],
          ", psit = ", solutions$psit[i],
          ", dist = ", solutions$dist[i]
        )
      ) +
      theme_bw()
  )

  dev.off()
}

# General helper functions used by pipeline R entrypoints.

# Call X/Y chromosome bins given autosomal purity and ploidy estimates.
callXchromsome = function(sex_lrrs, psi, psit, purity) {
  hap_adj_lrrs = log(((2^sex_lrrs) * psi) / (psi / 2), base = 2)
  hap_e = 1
  hap_psit = psit / 2
  hap_psi = (hap_e * (1 - purity)) + (purity * hap_psit)
  n = ((hap_psi * (2^hap_adj_lrrs)) - (hap_e * (1 - purity))) / purity

  return(n)
}

# Colour palette for integer copy-number states in ASCAT low-pass plots.
cols = c(
  c(
    "0" = "#1981be",
    "1" = "#56B4E9",
    "2" = "grey",
    "3" = "#E69F00",
    "4" = "#ffc342",
    "5" = "#FFAA42",
    "6" = "#FF9142",
    "7" = "#FF7742",
    "8" = "#FF5E42"
  ),
  rep("#FF4542", times = 100 - 8)
)
names(cols)[(9:100) + 1] = 9:100

# Expand segment rows to fixed-size genomic bins for MEDICC2 input preparation.
explode_ranges <- function(df, step_size) {
  exploded_list <- vector("list", nrow(df))

  for (i in 1:nrow(df)) {
    seq_starts <- seq(df$Start[i], df$End[i], by = step_size)
    seq_ends <- pmin(seq_starts + step_size - 1, df$End[i])

    exploded_list[[i]] <- data.frame(
      Chromosome = rep(df$Chromosome[i], length(seq_starts)),
      Start = seq_starts,
      End = seq_ends,
      Copies = rep(df$Copies[i], length(seq_starts))
    )
  }

  do.call(rbind, exploded_list)
}

# a function to return list of alternate solutions
make_solution <- function(res, purity, psit, label) {
    rho <- purity
    psi <- (2 * (1 - rho)) + (rho * psit)

    cn_auto <- ((psi * (2^segs_auto)) - (2 * (1 - rho))) / rho
    cn_auto_int <- round(cn_auto)
    cn_auto_int[cn_auto_int < 0] <- 0

    sex_segs <- expanded_segs[!autosome_index]
    cn_sex <- callXchromsome(
        sex_lrrs = sex_segs,
        psi = psi,
        psit = psit,
        purity = rho
    )
    cn_sex_int <- round(cn_sex)
    cn_sex_int[cn_sex_int < 0] <- 0

    out <- res
    out$CN <- c(cn_auto_int, cn_sex_int)
    out$contCN <- c(cn_auto, cn_sex)
    out$Psi <- psi
    out$PsiT <- psit
    out$Purity <- rho
    out$segs <- expanded_segs
    out$bins <- patient_lrr[, 5]
    out$sample_name <- id
    out$solution_label <- label
    out
}

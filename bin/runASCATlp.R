# ASCAT-style fitting on log2ratio only, with:
#  - fixed bug: fits loop now returns numbers (fit_mat is real)
#  - removed abs() bug
#  - optional soft penalty on excessive CN=0 usage during fitting (lambda0)

# lrrs: numeric vector of log2ratios
# psit: tumour ploidy to fit (e.g. 2 for diploid, 3 for triploid, etc.)

# Normal cells contribute copy number 2.
# Tumour cells contribute average copy number psit.

runASCATlp <- function(lrrs, fix_ploidy = 2, interval = 0.01,
                       min_purity = 0.2, max_purity = 1,
                       max_lrr = Inf, no_fit_psit = 2,
                       pad_ploidy = 0,
                       preset = FALSE,
                       preset_purity,
                       preset_ploidy,
                       lambda0 = 0.00) {   # <- set e.g. 0.01 to mildly discourage many CN=0

  gamma <- 1

  # Distance measure for a given (rho, psit)
  fit_lrr <- function(lrrs, rho, psit, gamma = 1, lambda0 = 0) {

    psi <- (2 * (1 - rho)) + (rho * psit)
    # psi: the average total copy number in the mixed sample
    # continuous CN
    n <- ((psi * (2^(lrrs / gamma))) - (2 * (1 - rho))) / rho
    # This converts each observed log2 ratio back into an inferred tumour copy number.

    # explanation:
    # observed relative CN = 2^lrr
    # mixed absolute CN    = psi * 2^lrr
    # tumour CN            = (mixed CN - normal contribution) / purity

    int_n <- round(n)
    int_n[int_n < 0] <- 0

    # SSE-to-integer + soft penalty for fraction of CN==0 calls
    fit <- mean((n - int_n)^2, na.rm = TRUE) + lambda0 * mean(int_n == 0, na.rm = TRUE)
    fit
  }

  # remove extreme lrrs for fitting (optional)
  lrrs_for_fit <- lrrs[is.finite(lrrs) & (lrrs < max_lrr)]

  if (!preset) {

    rhos  <- seq(min_purity, max_purity + interval, by = interval)
    psits <- seq(fix_ploidy - pad_ploidy, fix_ploidy + pad_ploidy, by = interval)

    # --- FIXED: actually RETURN numeric fits ---
    fits <- lapply(rhos, function(r) {
      sapply(psits, function(p) {
        fit_lrr(lrrs_for_fit, rho = r, psit = p, gamma = gamma, lambda0 = lambda0)
      })
    })

    fit_mat <- do.call(rbind, fits)
    colnames(fit_mat) <- psits
    rownames(fit_mat) <- rhos

    find_local_minima <- function(mat) {
      res <- NULL
      for (c in 2:(ncol(mat) - 1)) {
        for (r in 2:(nrow(mat) - 1)) {
          test_mat <- mat[(r - 1):(r + 1), (c - 1):(c + 1)]
          m <- which(test_mat == min(test_mat))
          if (m[1] == 5 && length(m) == 1) {
            hit <- c(colnames(test_mat)[2], rownames(test_mat)[2], test_mat[2, 2])
            res <- rbind(res, hit)
          }
        }
      }
      if (is.null(res)) {
        return(data.frame(psit = numeric(0), purity = numeric(0), dist = numeric(0)))
      }
      res <- data.frame(psit = as.numeric(res[, 1]),
                        purity = as.numeric(res[, 2]),
                        dist = round(as.numeric(res[, 3]), 5),
                        stringsAsFactors = FALSE)
      res[order(res$dist), ]
    }

    find_local_minima_single_ploidy <- function(mat) {
      res <- NULL
      for (r in 2:(nrow(mat) - 1)) {
        test_vec <- mat[(r - 1):(r + 1), 1]
        m <- which(test_vec == min(test_vec))
        if (m[1] == 2 && length(m) == 1) {
          hit <- c(psits, names(test_vec)[2], test_vec[2])
          res <- rbind(res, hit)
        }
      }
      if (is.null(res)) {
        return(data.frame(psit = numeric(0), purity = numeric(0), dist = numeric(0)))
      }
      res <- data.frame(psit = as.numeric(res[, 1]),
                        purity = as.numeric(res[, 2]),
                        dist = round(as.numeric(res[, 3]), 5),
                        stringsAsFactors = FALSE)
      res[order(res$dist), ]
    }

    lms <- if (ncol(fit_mat) == 1) {
         find_local_minima_single_ploidy(fit_mat)
     } else {
        find_local_minima(fit_mat)
    }

if (nrow(lms) == 0) {
  best_fit <- data.frame(psit = no_fit_psit, purity = min_purity, dist = Inf)
  near_best <- best_fit
} else {

  # Find near-equivalent best solutions
  best_dist <- min(lms$dist, na.rm = TRUE)

  near_best <- lms[lms$dist <= best_dist + 1e-4, , drop = FALSE]

  # Add diagnostics for tie-breaking
  near_best$delta_dist <- near_best$dist - best_dist
  near_best$ploidy_dev <- abs(near_best$psit - fix_ploidy)

  # Tie-break:
  # 1. best distance
  # 2. closest to expected/fixed ploidy
  # 3. higher purity
  near_best <- near_best[order(
    near_best$delta_dist,
    near_best$ploidy_dev,
    -near_best$purity
  ), ]

  best_fit <- near_best[1, , drop = FALSE]
}

  } else {
    best_fit <- data.frame(psit = preset_ploidy, purity = preset_purity, dist = NA)
    fit_mat  <- NULL
    lms      <- NULL
    near_best <- best_fit
  }

  # continuous CN for ALL lrrs under best fit
  rho <- best_fit$purity
  psi <- (2 * (1 - rho)) + (rho * best_fit$psit)
  n   <- ((psi * (2^(lrrs / gamma))) - (2 * (1 - rho))) / rho

  # integer CN calls
  n_int <- round(n)
  n_int[n_int < 0] <- 0
  list(
    CN        = n_int,
    Purity    = rho,
    Psi       = psi,
    PsiT      = best_fit$psit,
    contCN    = n,
    fit       = fit_mat,
    allsol    = lms,
    near_best = near_best
  )
}


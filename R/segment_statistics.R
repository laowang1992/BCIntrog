#' Summarize donor segments and genome proportions
#'
#' Summarize segment lengths by state and chromosome. The recurrent parent
#' background recovery rate uses the genome-equivalent definition
#' `(L_R + 0.5 * L_H) / L_total`; the donor genome equivalent is
#' `(L_D + 0.5 * L_H) / L_total`. Both are reported alongside the simpler
#' homozygous proportions.
#'
#' @details Segments with `Viterbi = NA` are reported as `Uninformative` and
#'   counted in `L_U` and `unassigned_ratio`; they are not forced into D, R or H.
#'
#' @param object An object of class `donor_segments` returned by
#'   [segment_detection()].
#' @param ... Additional arguments (unused).
#'
#' @return An object of class `summary_donor_segments`, a list with elements
#'   `overall`, `by_state`, `by_chrom`, `formulas`, `params` and `parents`.
#'
#' @examples
#' vcf <- system.file("extdata", "Demo.vcf.gz", package = "BCIntrog")
#' dat <- read_vcf(vcf)
#' sites <- select_site(dat, "9172", "Y633", "WS100")
#' seg <- segment_detection(sites, min_length = 500000)
#' summary(seg)
#'
#' @export
summary.donor_segments <- function(object, ...) {
  if (!inherits(object, "donor_segments")) {
    stop("`object` must be a `donor_segments` object.", call. = FALSE)
  }
  seg <- object$segments
  chr <- object$chr_length
  sites <- object$sites
  labels <- object$labels

  if (nrow(seg) == 0L) {
    stop("The `donor_segments` object contains no segments.", call. = FALSE)
  }

  state_key <- ifelse(is.na(seg$Viterbi), "Uninformative", seg$Viterbi)

  if (all(c("CHROM", "segment_id", "Confidence") %in% names(sites))) {
    site_key <- paste(sites$CHROM, sites$segment_id, sep = "\r")
    seg_key <- paste(seg$CHROM, seg$segment_id, sep = "\r")
    segment_conf <- vapply(seq_len(nrow(seg)), function(i) {
      idx <- which(site_key == seg_key[i])
      if (length(idx) == 0L || all(is.na(sites$Confidence[idx]))) {
        NA_real_
      } else {
        mean(sites$Confidence[idx], na.rm = TRUE)
      }
    }, numeric(1))
  } else if ("Mean_Confidence" %in% names(seg)) {
    segment_conf <- seg$Mean_Confidence
  } else {
    segment_conf <- rep(NA_real_, nrow(seg))
  }

  lengths <- tapply(seg$Length_bp, state_key, sum)
  get_len <- function(k) {
    if (k %in% names(lengths)) as.numeric(lengths[[k]]) else 0
  }
  L_D <- get_len("D")
  L_R <- get_len("R")
  L_H <- get_len("H")
  L_U <- get_len("Uninformative")
  L_total <- sum(as.numeric(chr$Length))
  L_assigned <- L_D + L_R + L_H

  prop <- function(x, denom) if (denom > 0) x / denom else NA_real_

  overall <- tibble::tibble(
    L_total = L_total,
    L_assigned = L_assigned,
    L_D = L_D,
    L_R = L_R,
    L_H = L_H,
    L_U = L_U,
    donor_ratio = prop(L_D, L_total),
    recurrent_ratio = prop(L_R, L_total),
    heterozygous_ratio = prop(L_H, L_total),
    donor_derived_ratio = prop(L_D + L_H, L_total),
    donor_genome_equivalent = prop(L_D + 0.5 * L_H, L_total),
    recurrent_recovery_allelic = prop(L_R + 0.5 * L_H, L_total),
    recurrent_recovery_homozygous = prop(L_R, L_total),
    unassigned_ratio = prop(L_U, L_total),
    n_segments = nrow(seg),
    n_donor_segments = sum(state_key == "D"),
    n_recurrent_segments = sum(state_key == "R"),
    n_heterozygous_segments = sum(state_key == "H"),
    n_uninformative_segments = sum(state_key == "Uninformative"),
    n_sites = nrow(sites),
    n_informative_sites = sum(!is.na(sites$Genotype)),
    missing_rate = if (nrow(sites) > 0) {
      mean(is.na(sites$Genotype))
    } else {
      NA_real_
    }
  )

  state_levels <- c("D", "R", "H", "Uninformative")
  state_labels <- c(
    D = labels$donor,
    R = labels$recurrent,
    H = labels$heterozygous,
    Uninformative = "Uninformative"
  )
  by_state <- tibble::tibble(
    Viterbi = c("D", "R", "H", NA_character_),
    State_Key = state_levels,
    State = unname(state_labels[state_levels])
  )
  by_state$Length_bp <- vapply(
    by_state$State_Key,
    function(k) get_len(k),
    numeric(1)
  )
  by_state$Proportion <- vapply(
    by_state$Length_bp,
    prop,
    numeric(1),
    denom = L_total
  )
  by_state$N_segments <- vapply(
    by_state$State_Key,
    function(k) sum(state_key == k),
    integer(1)
  )
  by_state$SNP_number <- vapply(
    by_state$State_Key,
    function(k) sum(seg$SNP_number[state_key == k], na.rm = TRUE),
    numeric(1)
  )
  by_state$Observed_SNP_number <- vapply(
    by_state$State_Key,
    function(k) sum(seg$Observed_SNP_number[state_key == k], na.rm = TRUE),
    numeric(1)
  )
  by_state$Mean_Confidence <- vapply(
    by_state$State_Key,
    function(k) {
      idx <- which(state_key == k)
      if (length(idx) == 0L) {
        return(NA_real_)
      }
      w <- seg$SNP_number[idx]
      x <- segment_conf[idx]
      ok <- !is.na(x) & !is.na(w) & w > 0
      if (!any(ok)) {
        return(NA_real_)
      }
      stats::weighted.mean(x[ok], w[ok])
    },
    numeric(1)
  )
  by_state <- tibble::as_tibble(lapply(by_state, unname))
  by_state <- by_state[, c(
    "Viterbi", "State", "Length_bp", "Proportion", "N_segments",
    "SNP_number", "Observed_SNP_number", "Mean_Confidence"
  )]

  by_chrom_list <- lapply(seq_len(nrow(chr)), function(i) {
    ch <- chr$CHROM[i]
    L <- as.numeric(chr$Length[i])
    sub <- seg[seg$CHROM == ch, , drop = FALSE]
    sk <- ifelse(is.na(sub$Viterbi), "Uninformative", sub$Viterbi)
    len <- function(k) sum(sub$Length_bp[sk == k], na.rm = TRUE)
    l_d <- len("D")
    l_r <- len("R")
    l_h <- len("H")
    l_u <- len("Uninformative")
    tibble::tibble(
      CHROM = ch,
      Length_bp = L,
      Donor_bp = l_d,
      Recurrent_bp = l_r,
      Heterozygous_bp = l_h,
      Uninformative_bp = l_u,
      donor_genome_equivalent = prop(l_d + 0.5 * l_h, L),
      recurrent_recovery_allelic = prop(l_r + 0.5 * l_h, L),
      recurrent_recovery_homozygous = prop(l_r, L),
      unassigned_ratio = prop(l_u, L),
      N_segments = nrow(sub),
      N_sites = sum(sites$CHROM == ch)
    )
  })
  by_chrom <- dplyr::bind_rows(by_chrom_list)

  formulas <- c(
    "L_total = sum(chr_length)",
    "donor_genome_equivalent = (L_D + 0.5 * L_H) / L_total",
    "recurrent_recovery_allelic = (L_R + 0.5 * L_H) / L_total",
    "recurrent_recovery_homozygous = L_R / L_total",
    "donor_ratio = L_D / L_total",
    "recurrent_ratio = L_R / L_total",
    "heterozygous_ratio = L_H / L_total",
    "unassigned_ratio = L_U / L_total"
  )

  structure(
    list(
      overall = overall,
      by_state = by_state,
      by_chrom = by_chrom,
      formulas = formulas,
      params = object$params,
      parents = labels
    ),
    class = "summary_donor_segments"
  )
}

#' @export
print.summary_donor_segments <- function(x, ...) {
  o <- x$overall
  pct <- function(z) ifelse(is.na(z), NA_character_, sprintf("%.2f%%", 100 * z))
  cat("BCIntrog segment summary\n")
  cat("  Total length:                 ",
      format(round(o$L_total / 1e6, 3), nsmall = 3), " Mb\n", sep = "")
  cat("  Donor genome equivalent:      ",
      pct(o$donor_genome_equivalent), "\n", sep = "")
  cat("  Recurrent recovery (allelic): ",
      pct(o$recurrent_recovery_allelic), "\n", sep = "")
  cat("  Recurrent homozygous:         ",
      pct(o$recurrent_recovery_homozygous), "\n", sep = "")
  cat("  Uninformative:                ",
      pct(o$unassigned_ratio), "\n", sep = "")
  cat("  Segments:                     ", o$n_segments, "\n", sep = "")
  cat("\nBy state:\n")
  print(as.data.frame(x$by_state), row.names = FALSE)
  invisible(x)
}

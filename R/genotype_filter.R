#' Select informative parental markers and call offspring genotypes
#'
#' Filter biallelic markers for which the donor and recurrent parents are
#' homozygous for opposite alleles, apply a parent GQ threshold, and recode the
#' offspring genotype as D (donor homozygous), R (recurrent homozygous), H
#' (heterozygous) or `NA` (missing/uninformative).
#'
#' @param data An object returned by [read_vcf()].
#' @param donor,recurrent,offspring Sample names for the donor parent, recurrent
#'   parent and backcross offspring.
#' @param GQ Parent GQ threshold. Only markers with both parent GQ values
#'   strictly greater than `GQ` are retained. Use `NULL` to skip GQ filtering.
#' @param biallelic_only Logical; retained for interface compatibility.
#'   Multiallelic records are not supported and are always removed.
#' @param verbose Logical; print progress messages.
#'
#' @return A tibble of class `donor_sites` with columns `CHROM`, `POS`, `REF`,
#'   `ALT`, `<donor>.GT`, `<recurrent>.GT`, `<offspring>.GT`, `Genotype`,
#'   `<donor>.GQ`, `<recurrent>.GQ` and `<offspring>.GQ`. The attributes
#'   `donor`, `recurrent`, `offspring`, `gq_threshold`, `has_gq`, `chrom_order`
#'   and `n_summary` are set for downstream functions.
#'
#' @examples
#' vcf <- system.file("extdata", "Demo.vcf.gz", package = "BCIntrog")
#' dat <- read_vcf(vcf)
#' sites <- select_site(dat, "9172", "Y633", "WS100")
#' sites
#'
#' @export
select_site <- function(data, donor, recurrent, offspring, GQ = 20,
                        biallelic_only = TRUE, verbose = TRUE) {
  if (!inherits(data, "donor_vcf")) {
    stop("`data` must be an object returned by `read_vcf()`.", call. = FALSE)
  }
  samples <- as.character(data$samples)
  requested <- c(donor, recurrent, offspring)
  if (anyNA(requested) || any(requested == "")) {
    stop("`donor`, `recurrent` and `offspring` must be non-empty sample names.",
         call. = FALSE)
  }
  if (anyDuplicated(requested)) {
    stop("`donor`, `recurrent` and `offspring` must be different samples.",
         call. = FALSE)
  }
  missing_samples <- setdiff(requested, samples)
  if (length(missing_samples) > 0L) {
    stop(
      "Samples not found in `data`: ",
      paste(missing_samples, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.null(GQ)) {
    .check_scalar_numeric(GQ, "GQ", min = 0)
  }

  sites <- data$sites
  gt <- data$gt
  gq <- data$gq
  n_total <- nrow(sites)

  biallelic <- !grepl(",", sites$ALT, fixed = TRUE)
  if (!all(biallelic)) {
    n_multi <- sum(!biallelic)
    if (verbose) {
      message("Removed ", n_multi, " multiallelic record(s).")
    }
    sites <- sites[biallelic, , drop = FALSE]
    gt <- gt[biallelic, , drop = FALSE]
    if (!is.null(gq)) {
      gq <- gq[biallelic, , drop = FALSE]
    }
  }
  n_biallelic <- nrow(sites)
  if (n_biallelic == 0L) {
    stop("No biallelic records remain.", call. = FALSE)
  }

  donor_p <- .parse_gt(gt[, donor])
  recurrent_p <- .parse_gt(gt[, recurrent])
  offspring_p <- .parse_gt(gt[, offspring])

  donor_dose <- donor_p$dosage
  recurrent_dose <- recurrent_p$dosage
  offspring_dose <- offspring_p$dosage

  parent_hom_opposite <-
    !is.na(donor_dose) & !is.na(recurrent_dose) &
    donor_dose %in% c(0L, 2L) &
    recurrent_dose %in% c(0L, 2L) &
    donor_dose != recurrent_dose

  gq_ok <- rep(TRUE, n_biallelic)
  gq_used <- FALSE
  n_gq_missing <- 0L
  if (!is.null(GQ)) {
    if (is.null(gq)) {
      warning(
        "GQ values are not available; GQ filtering will be skipped.",
        call. = FALSE
      )
    } else {
      gd <- gq[, donor]
      gr <- gq[, recurrent]
      donor_all_missing <- all(is.na(gd))
      recurrent_all_missing <- all(is.na(gr))
      if (donor_all_missing) {
        warning(
          "GQ values are entirely missing for the donor parent; ",
          "donor GQ filtering will be skipped.",
          call. = FALSE
        )
      }
      if (recurrent_all_missing) {
        warning(
          "GQ values are entirely missing for the recurrent parent; ",
          "recurrent GQ filtering will be skipped.",
          call. = FALSE
        )
      }
      if (donor_all_missing && recurrent_all_missing) {
        gq_ok <- rep(TRUE, n_biallelic)
      } else {
        gq_used <- TRUE
        gd_ok <- if (donor_all_missing) {
          rep(TRUE, n_biallelic)
        } else {
          !is.na(gd) & gd > GQ
        }
        gr_ok <- if (recurrent_all_missing) {
          rep(TRUE, n_biallelic)
        } else {
          !is.na(gr) & gr > GQ
        }
        gq_ok <- gd_ok & gr_ok
        n_gq_missing <- sum(
          parent_hom_opposite &
            ((is.na(gd) & !donor_all_missing) |
               (is.na(gr) & !recurrent_all_missing))
        )
        if (n_gq_missing > 0L) {
          warning(
            n_gq_missing,
            " informative marker(s) were removed because parent GQ was missing.",
            call. = FALSE
          )
        }
      }
    }
  }

  keep <- parent_hom_opposite & gq_ok

  # 后代基因型编码：与供体一致=D，与轮回亲本一致=R，杂合=H
  genotype <- rep(NA_character_, n_biallelic)
  d_idx <- which(
    !is.na(offspring_dose) & !is.na(donor_dose) &
      offspring_dose == donor_dose
  )
  r_idx <- which(
    !is.na(offspring_dose) & !is.na(recurrent_dose) &
      offspring_dose == recurrent_dose
  )
  genotype[d_idx] <- "D"
  genotype[r_idx] <- "R"
  h_idx <- which(
    !is.na(offspring_dose) & offspring_dose == 1L & is.na(genotype)
  )
  genotype[h_idx] <- "H"

  out <- tibble::tibble(
    CHROM = sites$CHROM[keep],
    POS = sites$POS[keep],
    REF = sites$REF[keep],
    ALT = sites$ALT[keep]
  )
  out[[.bt_gt_col(donor)]] <- gt[keep, donor]
  out[[.bt_gt_col(recurrent)]] <- gt[keep, recurrent]
  out[[.bt_gt_col(offspring)]] <- gt[keep, offspring]
  out$Genotype <- genotype[keep]
  if (!is.null(gq)) {
    out[[.bt_gq_col(donor)]] <- gq[keep, donor]
    out[[.bt_gq_col(recurrent)]] <- gq[keep, recurrent]
    out[[.bt_gq_col(offspring)]] <- gq[keep, offspring]
  }

  chrom_order <- data$chrom_order
  if (is.null(chrom_order)) {
    chrom_order <- unique(sites$CHROM)
  }
  ord <- order(match(out$CHROM, chrom_order), out$POS)
  out <- out[ord, , drop = FALSE]

  if (nrow(out) == 0L) {
    stop(
      "No informative markers remain after parental genotype and GQ filtering.",
      call. = FALSE
    )
  }

  attr(out, "donor") <- donor
  attr(out, "recurrent") <- recurrent
  attr(out, "offspring") <- offspring
  attr(out, "gq_threshold") <- GQ
  attr(out, "has_gq") <- gq_used
  attr(out, "chrom_order") <- chrom_order
  attr(out, "n_summary") <- list(
    n_total = n_total,
    n_biallelic = n_biallelic,
    n_parent_opposite = sum(parent_hom_opposite),
    n_gq_pass = sum(parent_hom_opposite & gq_ok),
    n_informative = nrow(out),
    n_offspring_missing = sum(is.na(out$Genotype)),
    n_invalid_gt = donor_p$n_invalid + recurrent_p$n_invalid +
      offspring_p$n_invalid,
    n_gq_missing = n_gq_missing
  )
  class(out) <- c("donor_sites", class(out))
  out
}

#' @export
print.donor_sites <- function(x, ...) {
  n <- attr(x, "n_summary")
  cat("BCIntrog informative markers\n")
  if (!is.null(n)) {
    cat("  Input records:       ", n$n_total, "\n", sep = "")
    cat("  Biallelic:           ", n$n_biallelic, "\n", sep = "")
    cat("  Opposite homozygous: ", n$n_parent_opposite, "\n", sep = "")
    cat("  Informative:         ", n$n_informative, "\n", sep = "")
    cat("  Offspring missing:   ", n$n_offspring_missing, "\n", sep = "")
  }
  cat("  Markers:             ", nrow(x), "\n", sep = "")
  invisible(x)
}

#' Detect donor-derived segments with a three-state HMM
#'
#' Identify genomic segments from HMM-inferred genetic states. Segments are
#' reported as 1-based closed intervals `[Start, End]`; the first segment starts
#' at 1, the last segment ends at the chromosome length, and adjacent segments
#' do not overlap (`Start[i + 1] = End[i] + 1`). Segments shorter than
#' `min_length` are merged or re-evaluated using the rules in `ALGORITHM.md`.
#'
#' @details If a chromosome has no observed offspring genotypes, or appears
#'   only in `chr_length`, a single `[1, L]` segment with `Viterbi = NA` and
#'   `State = "Uninformative"` is returned. Such segments are not forced into
#'   D, R or H.
#'
#' @param data A `donor_sites` object returned by [select_site()], or a data
#'   frame with columns `CHROM`, `POS` and `Genotype`.
#' @param min_length Minimum segment length in base pairs. Defaults to 500,000.
#'   Use 0 to skip small-segment processing.
#' @param chr_length Optional chromosome length table. Must have two columns:
#'   chromosome (character) and length (numeric). If `NULL`, the maximum marker
#'   position of each chromosome is used.
#' @param cross_length Expected distance between state changes, in base pairs.
#'   Defaults to 10,000,000.
#' @param emission A 3 x 3 emission probability matrix; each row must sum to 1.
#'   Defaults to [emission_probability()].
#' @param init_prob Numeric vector of initial state probabilities. Defaults to a
#'   uniform distribution.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `donor_segments`, a list with elements `sites`,
#'   `segments`, `chr_length`, `params` and `labels`. `segments` contains
#'   `CHROM`, `segment_id`, `Start`, `End`, `Length_bp`, `SNP_number`,
#'   `Observed_SNP_number`, `Viterbi` and `State`.
#'
#' @examples
#' vcf <- system.file("extdata", "Demo.vcf.gz", package = "BCIntrog")
#' dat <- read_vcf(vcf)
#' sites <- select_site(dat, "9172", "Y633", "WS100")
#' seg <- segment_detection(sites, min_length = 500000)
#' seg
#'
#' @export
segment_detection <- function(data, min_length = 500000, chr_length = NULL,
                              cross_length = 1e7, emission = NULL,
                              init_prob = NULL, verbose = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame, usually returned by `select_site()`.",
         call. = FALSE)
  }
  required <- c("CHROM", "POS", "Genotype")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0L) {
    stop(
      "`data` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (nrow(data) == 0L) {
    stop("`data` contains no markers.", call. = FALSE)
  }
  .check_scalar_numeric(min_length, "min_length", min = 0)
  .check_scalar_numeric(cross_length, "cross_length", min = 0)
  if (cross_length == 0) {
    stop("`cross_length` must be greater than 0.", call. = FALSE)
  }
  emission <- .validate_emission(emission)
  init_prob <- .validate_init_prob(init_prob)

  sites_all <- tibble::tibble(
    CHROM = as.character(data$CHROM),
    POS = as.integer(data$POS),
    Genotype = .validate_genotype(data$Genotype)
  )
  if (anyNA(sites_all$CHROM) || any(sites_all$CHROM == "")) {
    stop("`data$CHROM` contains missing or empty values.", call. = FALSE)
  }
  if (anyNA(sites_all$POS) || any(sites_all$POS <= 0L)) {
    stop("`data$POS` contains missing or non-positive values.", call. = FALSE)
  }
  if (anyDuplicated(paste(sites_all$CHROM, sites_all$POS, sep = ":"))) {
    stop("`data` contains duplicated CHROM/POS markers.", call. = FALSE)
  }

  chrom_order <- attr(data, "chrom_order")
  if (is.null(chrom_order)) {
    chrom_order <- unique(sites_all$CHROM)
  } else {
    chrom_order <- as.character(chrom_order)
    chrom_order <- chrom_order[chrom_order %in% sites_all$CHROM]
    extra <- setdiff(unique(sites_all$CHROM), chrom_order)
    chrom_order <- c(chrom_order, extra)
  }
  ord <- order(match(sites_all$CHROM, chrom_order), sites_all$POS)
  sites_all <- sites_all[ord, , drop = FALSE]

  labels <- .segment_labels(data)
  chr_tab <- .validate_chr_length(chr_length, sites_all, chrom_order)

  segments_list <- vector("list", nrow(chr_tab))
  sites_out_list <- vector("list", nrow(chr_tab))
  uninformative_chroms <- character(0)

  for (ci in seq_len(nrow(chr_tab))) {
    chrom <- chr_tab$CHROM[ci]
    L <- as.integer(chr_tab$Length[ci])
    sites <- sites_all[sites_all$CHROM == chrom, , drop = FALSE]
    n_sites <- nrow(sites)

    if (n_sites == 0L) {
      segments_list[[ci]] <- tibble::tibble(
        CHROM = chrom,
        segment_id = 1L,
        Start = 1L,
        End = L,
        Length_bp = L,
        SNP_number = 0L,
        Observed_SNP_number = 0L,
        Viterbi = NA_character_,
        State = "Uninformative",
        Mean_Prob_D = NA_real_,
        Mean_Prob_R = NA_real_,
        Mean_Prob_H = NA_real_,
        Mean_Confidence = NA_real_,
        unresolved = FALSE
      )
      sites_out_list[[ci]] <- sites
      next
    }

    if (all(is.na(sites$Genotype))) {
      uninformative_chroms <- c(uninformative_chroms, chrom)
      sites$Viterbi_initial <- NA_character_
      sites$Viterbi <- NA_character_
      sites$State <- "Uninformative"
      sites$Prob_D <- NA_real_
      sites$Prob_R <- NA_real_
      sites$Prob_H <- NA_real_
      sites$Confidence <- NA_real_
      sites$segment_id <- 1L
      sites_out_list[[ci]] <- sites
      segments_list[[ci]] <- tibble::tibble(
        CHROM = chrom,
        segment_id = 1L,
        Start = 1L,
        End = L,
        Length_bp = L,
        SNP_number = n_sites,
        Observed_SNP_number = 0L,
        Viterbi = NA_character_,
        State = "Uninformative",
        Mean_Prob_D = NA_real_,
        Mean_Prob_R = NA_real_,
        Mean_Prob_H = NA_real_,
        Mean_Confidence = NA_real_,
        unresolved = FALSE
      )
      next
    }

    hmm <- run_hmm(
      pos = sites$POS,
      observations = sites$Genotype,
      emission = emission,
      init_prob = init_prob,
      cross_length = cross_length
    )
    sites$Viterbi_initial <- hmm$Viterbi
    sites$Viterbi <- hmm$Viterbi
    sites$Prob_D <- hmm$Prob_D
    sites$Prob_R <- hmm$Prob_R
    sites$Prob_H <- hmm$Prob_H
    sites$Confidence <- hmm$Confidence

    seg <- .initial_segments(sites, hmm, L, labels)
    if (min_length > 0) {
      seg <- .process_small_segments(
        seg = seg,
        sites = sites,
        L = L,
        min_length = min_length,
        emission = emission,
        init_prob = init_prob,
        cross_length = cross_length,
        labels = labels,
        max_iter = 2L,
        chrom = chrom
      )
    }
    seg <- .recompute_segment_stats(seg, sites, labels)
    sites <- .assign_sites_to_segments(sites, seg, labels)
    segments_list[[ci]] <- seg
    sites_out_list[[ci]] <- sites
    if (verbose) {
      message(
        "Processed ", chrom, ": ", nrow(seg), " segment(s), ",
        n_sites, " markers."
      )
    }
  }

  if (length(uninformative_chroms) > 0L) {
    warning(
      "Chromosomes with no observed genotypes were reported as Uninformative: ",
      paste(uninformative_chroms, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  segments <- dplyr::bind_rows(segments_list)
  sites_out <- dplyr::bind_rows(sites_out_list)
  segments$CHROM <- as.character(segments$CHROM)
  segments$Viterbi <- as.character(segments$Viterbi)
  segments$State <- as.character(segments$State)
  segments <- tibble::as_tibble(
    segments[, c(
      "CHROM", "segment_id", "Start", "End", "Length_bp", "SNP_number",
      "Observed_SNP_number", "Viterbi", "State"
    ), drop = FALSE]
  )

  structure(
    list(
      sites = sites_out,
      segments = segments,
      chr_length = chr_tab,
      params = list(
        min_length = min_length,
        cross_length = cross_length,
        emission = emission,
        init_prob = init_prob,
        max_iter = 2L
      ),
      labels = labels
    ),
    class = c("donor_segments", "list")
  )
}

#' @export
print.donor_segments <- function(x, ...) {
  cat("BCIntrog segment detection result\n")
  cat("  Chromosomes: ", length(unique(x$segments$CHROM)), "\n", sep = "")
  cat("  Segments:    ", nrow(x$segments), "\n", sep = "")
  cat("  Markers:     ", nrow(x$sites), "\n", sep = "")
  cat("  min_length:  ", x$params$min_length, "\n", sep = "")
  cat("  cross_length:", x$params$cross_length, "\n", sep = "")
  invisible(x)
}

.segment_labels <- function(data) {
  donor <- attr(data, "donor")
  recurrent <- attr(data, "recurrent")
  if (is.null(donor) || is.null(recurrent) ||
      length(donor) != 1L || length(recurrent) != 1L ||
      is.na(donor) || is.na(recurrent)) {
    warning(
      "Donor/recurrent sample names were not found in the input attributes; ",
      "generic labels will be used.",
      call. = FALSE
    )
    return(list(
      donor = "Donor",
      recurrent = "Recurrent",
      heterozygous = "Heterozygous"
    ))
  }
  list(
    donor = as.character(donor),
    recurrent = as.character(recurrent),
    heterozygous = "Heterozygous"
  )
}

.initial_segments <- function(sites, hmm, L, labels) {
  # 由 Viterbi 游程生成初始片段；相邻转换标记的中点作为 breakpoint
  states <- hmm$Viterbi
  r <- rle(states)
  run_ends <- cumsum(r$lengths)
  run_starts <- c(1L, utils::head(run_ends, -1L) + 1L)
  k <- length(r$lengths)
  starts <- integer(k)
  ends <- integer(k)
  for (i in seq_len(k)) {
    starts[i] <- if (i == 1L) {
      1L
    } else {
      as.integer(floor(
        (sites$POS[run_starts[i] - 1L] + sites$POS[run_starts[i]]) / 2
      )) + 1L
    }
    ends[i] <- if (i == k) {
      as.integer(L)
    } else {
      as.integer(floor(
        (sites$POS[run_ends[i]] + sites$POS[run_ends[i] + 1L]) / 2
      ))
    }
  }
  seg <- tibble::tibble(
    CHROM = sites$CHROM[1L],
    segment_id = seq_len(k),
    Start = starts,
    End = ends,
    Length_bp = ends - starts + 1L,
    Viterbi = as.character(r$values),
    unresolved = FALSE
  )
  .recompute_segment_stats(seg, sites, labels)
}

.safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

.recompute_segment_stats <- function(seg, sites, labels) {
  if (nrow(seg) == 0L) {
    return(seg)
  }
  seg$Length_bp <- seg$End - seg$Start + 1L
  seg$SNP_number <- 0L
  seg$Observed_SNP_number <- 0L
  seg$Mean_Prob_D <- NA_real_
  seg$Mean_Prob_R <- NA_real_
  seg$Mean_Prob_H <- NA_real_
  seg$Mean_Confidence <- NA_real_
  for (i in seq_len(nrow(seg))) {
    idx <- which(sites$POS >= seg$Start[i] & sites$POS <= seg$End[i])
    seg$SNP_number[i] <- length(idx)
    if (length(idx) > 0L) {
      seg$Observed_SNP_number[i] <- sum(!is.na(sites$Genotype[idx]))
      seg$Mean_Prob_D[i] <- .safe_mean(sites$Prob_D[idx])
      seg$Mean_Prob_R[i] <- .safe_mean(sites$Prob_R[idx])
      seg$Mean_Prob_H[i] <- .safe_mean(sites$Prob_H[idx])
      seg$Mean_Confidence[i] <- .safe_mean(sites$Confidence[idx])
    }
  }
  seg$State <- .state_label(seg$Viterbi, labels)
  seg$segment_id <- seq_len(nrow(seg))
  seg <- seg[, c(
    "CHROM", "segment_id", "Start", "End", "Length_bp", "SNP_number",
    "Observed_SNP_number", "Viterbi", "State", "Mean_Prob_D", "Mean_Prob_R",
    "Mean_Prob_H", "Mean_Confidence", "unresolved"
  )]
  tibble::as_tibble(seg)
}

.assign_sites_to_segments <- function(sites, seg, labels) {
  sites$segment_id <- NA_integer_
  sites$Viterbi <- NA_character_
  sites$State <- "Uninformative"
  for (i in seq_len(nrow(seg))) {
    idx <- which(sites$POS >= seg$Start[i] & sites$POS <= seg$End[i])
    if (length(idx) > 0L) {
      sites$segment_id[idx] <- seg$segment_id[i]
      sites$Viterbi[idx] <- seg$Viterbi[i]
      sites$State[idx] <- seg$State[i]
    }
  }
  sites
}

.find_small_groups <- function(seg, min_length) {
  small <- which(
    seg$Length_bp < min_length & !is.na(seg$Viterbi) & !seg$unresolved
  )
  if (length(small) == 0L) {
    return(list())
  }
  breaks <- c(0L, which(diff(small) != 1L), length(small))
  groups <- vector("list", length(breaks) - 1L)
  for (i in seq_along(groups)) {
    idx <- small[(breaks[i] + 1L):breaks[i + 1L]]
    groups[[i]] <- c(min(idx), max(idx))
  }
  groups
}

.merge_adjacent <- function(seg, sites, labels) {
  if (nrow(seg) <= 1L) {
    return(.recompute_segment_stats(seg, sites, labels))
  }
  n <- nrow(seg)
  keep <- c(TRUE, logical(n - 1L))
  for (i in 2:n) {
    same_state <- (is.na(seg$Viterbi[i - 1L]) && is.na(seg$Viterbi[i])) ||
      (!is.na(seg$Viterbi[i - 1L]) && !is.na(seg$Viterbi[i]) &&
         seg$Viterbi[i - 1L] == seg$Viterbi[i])
    contiguous <- seg$Start[i] == seg$End[i - 1L] + 1L
    keep[i] <- !(same_state && contiguous)
  }
  grp <- cumsum(keep)
  new <- tibble::tibble(
    CHROM = seg$CHROM[1L],
    Start = as.integer(tapply(seg$Start, grp, min)),
    End = as.integer(tapply(seg$End, grp, max)),
    Viterbi = as.character(tapply(
      seg$Viterbi, grp,
      function(z) z[1L]
    )),
    unresolved = as.logical(tapply(seg$unresolved, grp, any))
  )
  .recompute_segment_stats(new, sites, labels)
}

.replace_segment_rows <- function(seg, from, to, new_rows) {
  before <- if (from > 1L) seg[seq_len(from - 1L), , drop = FALSE] else NULL
  after <- if (to < nrow(seg)) seg[(to + 1L):nrow(seg), , drop = FALSE] else NULL
  dplyr::bind_rows(before, new_rows, after)
}

.process_small_segments <- function(seg, sites, L, min_length, emission,
                                    init_prob, cross_length, labels,
                                    max_iter, chrom) {
  # 小片段处理：最多迭代 max_iter 轮，每轮先分组，再合并、切分或重跑 HMM
  if (nrow(seg) == 1L && seg$Length_bp[1L] < min_length) {
    return(seg)
  }
  for (iter in seq_len(max_iter)) {
    groups <- .find_small_groups(seg, min_length)
    if (length(groups) == 0L) {
      return(.recompute_segment_stats(seg, sites, labels))
    }
    sig_before <- paste(seg$Start, seg$End, seg$Viterbi, collapse = ";")
    n_small_before <- sum(
      seg$Length_bp < min_length & !is.na(seg$Viterbi) & !seg$unresolved
    )
    for (g in rev(groups)) {
      seg <- .process_small_group(
        seg = seg,
        g_start = g[1L],
        g_end = g[2L],
        sites = sites,
        L = L,
        min_length = min_length,
        emission = emission,
        init_prob = init_prob,
        cross_length = cross_length,
        labels = labels,
        chrom = chrom
      )
    }
    seg <- .merge_adjacent(seg, sites, labels)
    sig_after <- paste(seg$Start, seg$End, seg$Viterbi, collapse = ";")
    n_small_after <- sum(
      seg$Length_bp < min_length & !is.na(seg$Viterbi) & !seg$unresolved
    )
    if (n_small_after >= n_small_before && identical(sig_before, sig_after)) {
      break
    }
  }
  small_left <- seg$Length_bp < min_length & !is.na(seg$Viterbi) &
    !seg$unresolved
  if (any(small_left)) {
    seg$unresolved[small_left] <- TRUE
    warning(
      "Chromosome ",
      chrom,
      ": some segments remain shorter than `min_length` after ",
      max_iter,
      " iteration(s). They are retained in the output.",
      call. = FALSE
    )
  }
  .recompute_segment_stats(seg, sites, labels)
}

.process_small_group <- function(seg, g_start, g_end, sites, L, min_length,
                                 emission, init_prob, cross_length, labels,
                                 chrom) {
  group <- seg[g_start:g_end, , drop = FALSE]
  left <- if (g_start > 1L) seg[g_start - 1L, , drop = FALSE] else NULL
  right <- if (g_end < nrow(seg)) seg[g_end + 1L, , drop = FALSE] else NULL
  W <- group$End[nrow(group)] - group$Start[1L] + 1L

  if (is.null(left) && is.null(right)) {
    if (nrow(group) <= 1L) {
      return(seg)
    }
    return(.rerun_group_hmm(
      seg, g_start, g_end, group, sites, emission, init_prob, cross_length,
      labels
    ))
  }

  if (W >= min_length) {
    return(.rerun_group_hmm(
      seg, g_start, g_end, group, sites, emission, init_prob, cross_length,
      labels
    ))
  }

  if (!is.null(left) && !is.null(right)) {
    if (identical(left$Viterbi, right$Viterbi)) {
      new_row <- tibble::tibble(
        CHROM = chrom,
        Start = left$Start,
        End = right$End,
        Viterbi = left$Viterbi,
        unresolved = FALSE
      )
      return(.replace_segment_rows(seg, g_start - 1L, g_end + 1L, new_row))
    }
    split <- .split_group_bounds(group, sites, W)
    if (is.null(split)) {
      warning(
        "Chromosome ",
        chrom,
        ": a short segment could not be split by SNP count; it was merged ",
        "into the left flank.",
        call. = FALSE
      )
      new_row <- tibble::tibble(
        CHROM = chrom,
        Start = left$Start,
        End = group$End[nrow(group)],
        Viterbi = left$Viterbi,
        unresolved = FALSE
      )
      return(.replace_segment_rows(seg, g_start - 1L, g_end, new_row))
    }
    new_rows <- tibble::tibble(
      CHROM = chrom,
      Start = c(left$Start, split + 1L),
      End = c(split, right$End),
      Viterbi = c(left$Viterbi, right$Viterbi),
      unresolved = FALSE
    )
    return(.replace_segment_rows(seg, g_start - 1L, g_end + 1L, new_rows))
  }

  if (!is.null(left)) {
    new_row <- tibble::tibble(
      CHROM = chrom,
      Start = left$Start,
      End = group$End[nrow(group)],
      Viterbi = left$Viterbi,
      unresolved = FALSE
    )
    return(.replace_segment_rows(seg, g_start - 1L, g_end, new_row))
  }

  new_row <- tibble::tibble(
    CHROM = chrom,
    Start = group$Start[1L],
    End = right$End,
    Viterbi = right$Viterbi,
    unresolved = FALSE
  )
  .replace_segment_rows(seg, g_start, g_end + 1L, new_row)
}

.split_group_bounds <- function(group, sites, W) {
  left_len <- floor(W / 2)
  if (left_len > 0L && left_len < W) {
    return(as.integer(group$Start[1L] + left_len - 1L))
  }
  p <- sites$POS[
    sites$POS >= group$Start[1L] & sites$POS <= group$End[nrow(group)]
  ]
  if (length(p) < 2L) {
    return(NULL)
  }
  k <- floor(length(p) / 2)
  split <- as.integer(floor((p[k] + p[k + 1L]) / 2))
  if (split < group$Start[1L] || split >= group$End[nrow(group)]) {
    return(NULL)
  }
  split
}

.rerun_group_hmm <- function(seg, g_start, g_end, group, sites, emission,
                             init_prob, cross_length, labels) {
  # 将小片段视为 SNP：位置取片段中点，观测值取上一轮 HMM 的 Viterbi 状态
  marker_pos <- (group$Start + group$End) / 2
  marker_obs <- group$Viterbi
  marker_obs[group$Observed_SNP_number == 0L] <- NA_character_
  if (all(is.na(marker_obs))) {
    warning(
      "A short-segment group has no observed markers; it was left unchanged.",
      call. = FALSE
    )
    return(seg)
  }
  hmm <- run_hmm(
    pos = marker_pos,
    observations = marker_obs,
    emission = emission,
    init_prob = init_prob,
    cross_length = cross_length
  )
  new_state <- hmm$Viterbi
  new_state[is.na(new_state)] <- group$Viterbi[is.na(new_state)]
  group$Viterbi <- new_state
  group$unresolved <- FALSE
  .replace_segment_rows(seg, g_start, g_end, group)
}

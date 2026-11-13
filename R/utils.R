# Internal helpers -----------------------------------------------------------

.bt_states <- c("D", "R", "H")

.bt_gt_col <- function(sample) paste0(sample, ".GT")

.bt_gq_col <- function(sample) paste0(sample, ".GQ")

.check_scalar_numeric <- function(x, name, min = NULL, max = NULL,
                                  finite = TRUE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be a single non-missing numeric value.",
         call. = FALSE)
  }
  if (finite && !is.finite(x)) {
    stop("`", name, "` must be finite.", call. = FALSE)
  }
  if (!is.null(min) && x < min) {
    stop("`", name, "` must be >= ", min, ".", call. = FALSE)
  }
  if (!is.null(max) && x > max) {
    stop("`", name, "` must be <= ", max, ".", call. = FALSE)
  }
  if (integer && x != floor(x)) {
    stop("`", name, "` must be a whole number.", call. = FALSE)
  }
  invisible(x)
}

.log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

.parse_gt <- function(x) {
  x <- chartr("|", "/", as.character(x))
  out <- rep(NA_integer_, length(x))
  out[x == "0/0"] <- 0L
  out[x %in% c("0/1", "1/0")] <- 1L
  out[x == "1/1"] <- 2L
  missing <- is.na(x) | x %in% c("./.", ".", "")
  list(
    dosage = out,
    n_missing = sum(missing),
    n_invalid = sum(!missing & is.na(out))
  )
}

.validate_genotype <- function(x, name = "Genotype") {
  if (!is.character(x) && !is.factor(x)) {
    stop("`", name, "` must be a character or factor vector.", call. = FALSE)
  }
  x <- as.character(x)
  allowed <- c(.bt_states, NA_character_)
  bad <- !(x %in% allowed)
  if (any(bad)) {
    stop(
      "`", name, "` contains invalid values: ",
      paste(unique(x[bad]), collapse = ", "),
      ". Allowed values are D, R, H or NA.",
      call. = FALSE
    )
  }
  x
}

.validate_emission <- function(emission = NULL) {
  if (is.null(emission)) {
    return(.default_emission())
  }
  if (is.data.frame(emission)) {
    emission <- as.matrix(emission)
  }
  if (!is.matrix(emission) || !is.numeric(emission) ||
      nrow(emission) != 3L || ncol(emission) != 3L) {
    stop("`emission` must be a 3 x 3 numeric matrix.", call. = FALSE)
  }
  if (anyNA(emission) || any(!is.finite(emission))) {
    stop("`emission` must not contain NA, NaN or infinite values.",
         call. = FALSE)
  }
  if (any(emission < 0)) {
    stop("`emission` must contain non-negative probabilities.", call. = FALSE)
  }
  rn <- rownames(emission)
  cn <- colnames(emission)
  if (!is.null(rn) || !is.null(cn)) {
    if (is.null(rn) || is.null(cn) ||
        !all(.bt_states %in% rn) || !all(.bt_states %in% cn)) {
      stop("`emission` dimnames must contain D, R and H for both rows and columns.",
           call. = FALSE)
    }
    emission <- emission[.bt_states, .bt_states, drop = FALSE]
  } else {
    dimnames(emission) <- list(.bt_states, .bt_states)
  }
  rs <- rowSums(emission)
  if (any(abs(rs - 1) > 1e-8)) {
    stop("Each row of `emission` must sum to 1.", call. = FALSE)
  }
  emission
}

.validate_init_prob <- function(init_prob = NULL) {
  if (is.null(init_prob)) {
    init_prob <- stats::setNames(rep(1 / 3, 3), .bt_states)
  }
  if (!is.numeric(init_prob) || length(init_prob) != 3L) {
    stop("`init_prob` must be a numeric vector of length 3.", call. = FALSE)
  }
  if (anyNA(init_prob) || any(!is.finite(init_prob)) || any(init_prob < 0)) {
    stop("`init_prob` must contain finite non-negative values.", call. = FALSE)
  }
  nms <- names(init_prob)
  if (!is.null(nms)) {
    if (!all(.bt_states %in% nms)) {
      stop("`init_prob` names must include D, R and H.", call. = FALSE)
    }
    init_prob <- init_prob[.bt_states]
  } else {
    names(init_prob) <- .bt_states
  }
  if (abs(sum(init_prob) - 1) > 1e-8) {
    stop("`init_prob` must sum to 1.", call. = FALSE)
  }
  init_prob
}

.validate_chr_length <- function(chr_length, sites, chrom_order) {
  pos_max <- tapply(sites$POS, sites$CHROM, max)
  default <- tibble::tibble(
    CHROM = chrom_order,
    Length = as.numeric(pos_max[chrom_order])
  )

  if (is.null(chr_length)) {
    return(default)
  }

  if (is.data.frame(chr_length)) {
    tab <- as.data.frame(chr_length, stringsAsFactors = FALSE)
  } else if (is.matrix(chr_length)) {
    tab <- as.data.frame(chr_length, stringsAsFactors = FALSE)
  } else {
    stop("`chr_length` must be NULL, a data.frame, a tibble or a matrix.",
         call. = FALSE)
  }
  if (ncol(tab) != 2L) {
    stop("`chr_length` must have exactly two columns: chromosome and length.",
         call. = FALSE)
  }
  chroms <- tab[[1L]]
  lengths <- tab[[2L]]
  if (!is.character(chroms) && !is.factor(chroms)) {
    stop("The first column of `chr_length` must be character.",
         call. = FALSE)
  }
  chroms <- as.character(chroms)
  if (!is.numeric(lengths)) {
    stop("The second column of `chr_length` must be numeric.", call. = FALSE)
  }
  if (anyNA(chroms) || any(chroms == "")) {
    stop("`chr_length` contains missing or empty chromosome names.",
         call. = FALSE)
  }
  if (anyDuplicated(chroms)) {
    stop("`chr_length` contains duplicated chromosome names.", call. = FALSE)
  }
  if (anyNA(lengths) || any(!is.finite(lengths)) || any(lengths <= 0)) {
    stop("`chr_length` lengths must be finite positive numbers.",
         call. = FALSE)
  }
  if (any(lengths != floor(lengths))) {
    stop("`chr_length` lengths must be whole numbers.", call. = FALSE)
  }

  user <- tibble::tibble(CHROM = chroms, Length = as.numeric(lengths))
  missing_chrom <- setdiff(chrom_order, user$CHROM)
  if (length(missing_chrom) > 0L) {
    warning(
      "Some chromosomes are missing from `chr_length`; using the maximum SNP position: ",
      paste(missing_chrom, collapse = ", "),
      call. = FALSE
    )
    user <- dplyr::bind_rows(user, default[default$CHROM %in% missing_chrom, ])
  }

  bad <- user$CHROM %in% chrom_order &
    user$Length < as.numeric(pos_max[user$CHROM])
  if (any(bad)) {
    stop(
      "`chr_length` is shorter than the maximum SNP position for: ",
      paste(user$CHROM[bad], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  extra <- setdiff(user$CHROM, chrom_order)
  ordered <- user[user$CHROM %in% chrom_order, ]
  ordered <- ordered[match(chrom_order, ordered$CHROM), ]
  if (length(extra) > 0L) {
    ordered <- dplyr::bind_rows(ordered, user[user$CHROM %in% extra, ])
  }
  tibble::as_tibble(ordered)
}

.state_label <- function(viterbi, labels) {
  out <- rep(NA_character_, length(viterbi))
  out[viterbi == "D"] <- labels$donor
  out[viterbi == "R"] <- labels$recurrent
  out[viterbi == "H"] <- labels$heterozygous
  out[is.na(viterbi)] <- "Uninformative"
  out
}

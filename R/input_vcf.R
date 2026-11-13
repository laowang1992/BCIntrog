#' Read a VCF or VCF.GZ file
#'
#' Read a VCF file with \pkg{vcfR} and extract site information, raw genotypes
#' (`GT`) and, when available, genotype quality (`GQ`) for the requested
#' samples. Only biallelic records are retained by default; both SNPs and
#' indels are allowed.
#'
#' @param file Path to a VCF or VCF.GZ file.
#' @param samples Character vector of sample names to retain. `NULL` keeps all
#'   samples.
#' @param fields FORMAT fields to extract. Must include `"GT"`; `"GQ"` is
#'   optional.
#' @param biallelic_only Logical; if `TRUE`, remove records whose `ALT` contains
#'   a comma (multiallelic records).
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `donor_vcf`, a list with elements:
#'   \describe{
#'     \item{sites}{A tibble with `CHROM`, `POS`, `REF` and `ALT`.}
#'     \item{gt}{A character matrix of raw genotypes (sites x samples).}
#'     \item{gq}{A numeric matrix of GQ values, or `NULL`.}
#'     \item{samples}{Sample names.}
#'     \item{has_gq}{Logical; whether usable GQ values are present.}
#'     \item{chrom_order}{Chromosome order as they first appear in the VCF.}
#'     \item{source_file}{The input file path.}
#'     \item{n_raw}{Number of records read before biallelic filtering.}
#'     \item{n_biallelic}{Number of records after biallelic filtering.}
#'   }
#'
#' @examples
#' vcf <- system.file("extdata", "Demo.vcf.gz", package = "BCIntrog")
#' dat <- read_vcf(vcf)
#' dat
#'
#' @export
read_vcf <- function(file, samples = NULL, fields = c("GT", "GQ"),
                     biallelic_only = TRUE, verbose = FALSE) {
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop("`file` must be a single file path.", call. = FALSE)
  }
  if (!file.exists(file) || dir.exists(file)) {
    stop("VCF file does not exist: ", file, call. = FALSE)
  }
  fields <- unique(as.character(fields))
  if (!"GT" %in% fields) {
    stop("`fields` must include \"GT\".", call. = FALSE)
  }
  unsupported <- setdiff(fields, c("GT", "GQ"))
  if (length(unsupported) > 0L) {
    stop(
      "Unsupported `fields`: ", paste(unsupported, collapse = ", "),
      ". Only GT and GQ are supported.",
      call. = FALSE
    )
  }
  want_gq <- "GQ" %in% fields

  vcf <- vcfR::read.vcfR(file, verbose = verbose)
  fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  if (nrow(fix) == 0L) {
    stop("The VCF contains no records.", call. = FALSE)
  }

  gt <- vcfR::extract.gt(vcf, element = "GT")
  if (is.null(dim(gt))) {
    gt <- matrix(gt, ncol = 1L)
    colnames(gt) <- setdiff(colnames(vcf@gt), "FORMAT")
  }
  gt <- as.matrix(gt)
  sample_names <- colnames(gt)
  if (is.null(sample_names) || length(sample_names) == 0L) {
    stop("The VCF contains no samples.", call. = FALSE)
  }
  if (anyDuplicated(sample_names)) {
    stop("The VCF contains duplicated sample names.", call. = FALSE)
  }

  if (is.null(samples)) {
    samples <- sample_names
  }
  samples <- as.character(samples)
  missing_samples <- setdiff(samples, sample_names)
  if (length(missing_samples) > 0L) {
    stop(
      "Samples not found in the VCF: ",
      paste(missing_samples, collapse = ", "),
      ". Available samples: ",
      paste(sample_names, collapse = ", "),
      call. = FALSE
    )
  }
  sample_idx <- match(samples, sample_names)
  gt <- gt[, sample_idx, drop = FALSE]

  has_gq_field <- FALSE
  if (want_gq) {
    format_strings <- as.character(vcf@gt[, "FORMAT"])
    format_fields <- strsplit(format_strings, ":", fixed = TRUE)
    has_gq_field <- any(vapply(
      format_fields,
      function(z) "GQ" %in% z,
      logical(1L)
    ))
  }

  gq <- NULL
  has_gq <- FALSE
  if (want_gq && has_gq_field) {
    gq <- suppressWarnings(
      vcfR::extract.gt(vcf, element = "GQ", as.numeric = TRUE)
    )
    if (!is.null(gq)) {
      if (is.null(dim(gq))) {
        gq <- matrix(gq, ncol = 1L)
        colnames(gq) <- sample_names
      }
      gq <- as.matrix(gq)[, sample_idx, drop = FALSE]
      has_gq <- any(!is.na(gq))
    }
    if (!has_gq) {
      warning(
        "The VCF contains a GQ field but all GQ values are missing; ",
        "GQ filtering will be skipped.",
        call. = FALSE
      )
      gq <- NULL
    }
  } else if (want_gq) {
    warning(
      "The VCF does not contain a GQ field; GQ filtering will be skipped.",
      call. = FALSE
    )
  }

  chrom <- as.character(fix$CHROM)
  pos <- suppressWarnings(as.integer(fix$POS))
  ref <- as.character(fix$REF)
  alt <- as.character(fix$ALT)
  n_raw <- length(chrom)

  if (anyNA(chrom) || any(chrom == "")) {
    stop("The VCF contains missing or empty CHROM values.", call. = FALSE)
  }
  if (anyNA(pos) || any(pos <= 0L)) {
    stop("The VCF contains missing, non-numeric or non-positive POS values.",
         call. = FALSE)
  }
  if (anyNA(ref) || anyNA(alt)) {
    stop("The VCF contains missing REF or ALT values.", call. = FALSE)
  }

  keep <- rep(TRUE, n_raw)
  if (isTRUE(biallelic_only)) {
    keep <- !grepl(",", alt, fixed = TRUE)
    n_drop <- sum(!keep)
    if (n_drop > 0L && verbose) {
      message("Removed ", n_drop, " multiallelic record(s).")
    }
  }
  chrom <- chrom[keep]
  pos <- pos[keep]
  ref <- ref[keep]
  alt <- alt[keep]
  gt <- gt[keep, , drop = FALSE]
  if (!is.null(gq)) {
    gq <- gq[keep, , drop = FALSE]
  }
  n_biallelic <- length(chrom)

  if (n_biallelic == 0L) {
    stop("No records remain after biallelic filtering.", call. = FALSE)
  }

  dup_key <- paste(chrom, pos, sep = ":")
  if (anyDuplicated(dup_key)) {
    dup <- unique(dup_key[duplicated(dup_key)])
    stop(
      "Duplicated CHROM/POS records are not allowed. Examples: ",
      paste(utils::head(dup, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  sites <- tibble::tibble(
    CHROM = chrom,
    POS = pos,
    REF = ref,
    ALT = alt
  )

  structure(
    list(
      sites = sites,
      gt = gt,
      gq = gq,
      samples = samples,
      has_gq = isTRUE(has_gq),
      chrom_order = unique(chrom),
      source_file = file,
      n_raw = n_raw,
      n_biallelic = n_biallelic
    ),
    class = c("donor_vcf", "list")
  )
}

#' @export
print.donor_vcf <- function(x, ...) {
  cat("BCIntrog VCF object\n")
  cat("  File:    ", x$source_file, "\n", sep = "")
  cat("  Records: ", x$n_biallelic, " biallelic (", x$n_raw,
      " raw)\n", sep = "")
  cat("  Samples: ", length(x$samples), "\n", sep = "")
  cat("  GQ:      ", if (x$has_gq) "available" else "not available", "\n",
      sep = "")
  invisible(x)
}

#' @export
as.data.frame.donor_vcf <- function(x, row.names = NULL, optional = FALSE,
                                    ...) {
  out <- as.data.frame(x$sites, stringsAsFactors = FALSE)
  for (s in x$samples) {
    out[[.bt_gt_col(s)]] <- x$gt[, s]
    if (!is.null(x$gq)) {
      out[[.bt_gq_col(s)]] <- x$gq[, s]
    }
  }
  out
}

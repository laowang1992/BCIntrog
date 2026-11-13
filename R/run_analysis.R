#' Run the complete BCIntrog analysis pipeline
#'
#' Read a VCF, select informative parental markers, call offspring genotypes,
#' detect donor segments, summarize the results and create a genome-wide plot.
#' This is a convenience wrapper around [read_vcf()], [select_site()],
#' [segment_detection()], [summary()] and [visualization()].
#'
#' @details Returned segments may contain `Uninformative` segments for
#'   chromosomes without observed offspring genotypes.
#'
#' @param file Path to a VCF or VCF.GZ file.
#' @param donor,recurrent,offspring Sample names for the donor parent, recurrent
#'   parent and backcross offspring.
#' @param GQ Parent GQ threshold. Use `NULL` to skip GQ filtering.
#' @param min_length Minimum segment length in base pairs.
#' @param chr_length Optional chromosome length table.
#' @param cross_length Expected distance between state changes, in base pairs.
#' @param emission Optional 3 x 3 emission probability matrix; each row must
#'   sum to 1.
#' @param init_prob Optional vector of initial state probabilities.
#' @param chr,label Optional chromosome selection and display labels for the
#'   plot.
#' @param output_dir Optional directory for result files. If `NULL`, no files
#'   are written.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `bcintrog_result`, a list with elements `vcf`,
#'   `sites`, `segments`, `summary` and `plot`.
#'
#' @examples
#' vcf <- system.file("extdata", "Demo.vcf.gz", package = "BCIntrog")
#' res <- run_analysis(vcf, "9172", "Y633", "WS100")
#' res$summary
#'
#' @export
run_analysis <- function(file, donor, recurrent, offspring, GQ = 20,
                         min_length = 500000, chr_length = NULL,
                         cross_length = 1e7, emission = NULL,
                         init_prob = NULL, chr = NULL, label = NULL,
                         output_dir = NULL, verbose = TRUE) {
  dat <- read_vcf(
    file = file,
    samples = unique(c(donor, recurrent, offspring)),
    verbose = verbose
  )
  sites <- select_site(
    data = dat,
    donor = donor,
    recurrent = recurrent,
    offspring = offspring,
    GQ = GQ,
    verbose = verbose
  )
  seg <- segment_detection(
    data = sites,
    min_length = min_length,
    chr_length = chr_length,
    cross_length = cross_length,
    emission = emission,
    init_prob = init_prob,
    verbose = verbose
  )
  summ <- summary(seg)
  p <- visualization(seg, chr = chr, label = label)

  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L ||
        is.na(output_dir)) {
      stop("`output_dir` must be a single directory path.", call. = FALSE)
    }
    if (!dir.exists(output_dir)) {
      ok <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if (!ok) {
        stop("Could not create `output_dir`: ", output_dir, call. = FALSE)
      }
    }
    utils::write.table(
      as.data.frame(sites),
      file = file.path(output_dir, "informative_snps.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )
    utils::write.table(
      as.data.frame(seg$sites),
      file = file.path(output_dir, "hmm_sites.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )
    utils::write.table(
      as.data.frame(seg$segments),
      file = file.path(output_dir, "segments.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )
    utils::write.table(
      as.data.frame(summ$overall),
      file = file.path(output_dir, "summary_overall.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )
    utils::write.table(
      as.data.frame(summ$by_state),
      file = file.path(output_dir, "summary_by_state.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )
    utils::write.table(
      as.data.frame(summ$by_chrom),
      file = file.path(output_dir, "summary_by_chrom.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )
    ggplot2::ggsave(
      filename = file.path(output_dir, "genome_segments.pdf"),
      plot = p, width = 12, height = 8
    )
    ggplot2::ggsave(
      filename = file.path(output_dir, "genome_segments.png"),
      plot = p, width = 12, height = 8, dpi = 300
    )
  }

  structure(
    list(
      vcf = dat,
      sites = sites,
      segments = seg,
      summary = summ,
      plot = p
    ),
    class = c("bcintrog_result", "list")
  )
}

#' @export
print.bcintrog_result <- function(x, ...) {
  cat("BCIntrog analysis result\n")
  cat("  Informative markers: ", nrow(x$sites), "\n", sep = "")
  cat("  Segments:            ", nrow(x$segments$segments), "\n", sep = "")
  cat("  Donor equivalent:    ",
      sprintf("%.2f%%", 100 * x$summary$overall$donor_genome_equivalent),
      "\n", sep = "")
  cat("  Recurrent recovery:  ",
      sprintf("%.2f%%", 100 * x$summary$overall$recurrent_recovery_allelic),
      "\n", sep = "")
  invisible(x)
}

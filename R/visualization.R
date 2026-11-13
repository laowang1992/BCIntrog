#' Plot genome-wide marker genotypes and inferred segments
#'
#' Draw a genome-wide plot with one horizontal bar per chromosome. Recurrent
#' parent segments are shown as the grey background; donor (D) and heterozygous
#' (H) segments are coloured. Marker genotypes are drawn as short vertical
#' ticks coloured by the observed genotype.
#'
#' @param data An object of class `donor_segments` returned by
#'   [segment_detection()].
#' @param chr Character vector of chromosomes to display, in display order.
#'   Defaults to all chromosomes in `data`.
#' @param label Character vector of chromosome labels. Defaults to `chr`. If
#'   longer than `chr`, extra labels are ignored with a warning; if shorter,
#'   missing labels are filled with the original chromosome names.
#' @param colors Named character vector of colours for `D`, `R` and `H`.
#' @param file Optional output file (`.pdf` or `.png`). If `NULL`, the plot is
#'   returned without saving.
#' @param width,height Plot dimensions in inches.
#' @param dpi Resolution for raster output and optional rasterisation.
#' @param rasterize Logical; if `TRUE` and \pkg{ggrastr} is installed, marker
#'   ticks are rasterised.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' vcf <- system.file("extdata", "Demo.vcf.gz", package = "BCIntrog")
#' dat <- read_vcf(vcf)
#' sites <- select_site(dat, "9172", "Y633", "WS100")
#' seg <- segment_detection(sites, min_length = 500000)
#' p <- visualization(seg)
#'
#' @export
visualization <- function(data, chr = NULL, label = NULL,
                          colors = c(D = "red", R = "blue", H = "green"),
                          file = NULL, width = 12, height = 8, dpi = 300,
                          rasterize = TRUE, ...) {
  if (!inherits(data, "donor_segments")) {
    stop("`data` must be a `donor_segments` object.", call. = FALSE)
  }
  chr_tab <- data$chr_length
  if (nrow(chr_tab) == 0L) {
    stop("`data` contains no chromosome lengths.", call. = FALSE)
  }
  all_chrom <- as.character(chr_tab$CHROM)

  if (is.null(chr)) {
    chr <- all_chrom
  } else {
    if (!is.character(chr) || length(chr) == 0L || anyNA(chr)) {
      stop("`chr` must be a non-empty character vector without NA.",
           call. = FALSE)
    }
    if (anyDuplicated(chr)) {
      stop("`chr` must not contain duplicated chromosome names.",
           call. = FALSE)
    }
    missing_chr <- setdiff(chr, all_chrom)
    if (length(missing_chr) > 0L) {
      stop(
        "Chromosomes not found in `data`: ",
        paste(missing_chr, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (is.null(label)) {
    label <- chr
  } else {
    if (!is.character(label)) {
      stop("`label` must be a character vector.", call. = FALSE)
    }
    if (length(label) > length(chr)) {
      warning(
        "`label` is longer than `chr`; extra labels were ignored.",
        call. = FALSE
      )
      label <- label[seq_along(chr)]
    } else if (length(label) < length(chr)) {
      warning(
        "`label` is shorter than `chr`; missing labels were filled with ",
        "the original chromosome names.",
        call. = FALSE
      )
      label <- c(label, chr[(length(label) + 1L):length(chr)])
    }
    if (anyNA(label)) {
      warning(
        "`label` contains NA; missing labels were filled with the original ",
        "chromosome names.",
        call. = FALSE
      )
      label[is.na(label)] <- chr[is.na(label)]
    }
  }

  if (is.null(names(colors))) {
    if (length(colors) != 3L) {
      stop("`colors` must be a named vector with entries D, R and H.",
           call. = FALSE)
    }
    names(colors) <- c("D", "R", "H")
  }
  if (!all(c("D", "R", "H") %in% names(colors))) {
    stop("`colors` must contain entries named D, R and H.", call. = FALSE)
  }
  colors <- colors[c("D", "R", "H")]

  chr_df <- tibble::tibble(
    CHROM = chr,
    Label = label,
    y = rev(seq_along(chr) * 3)
  )
  idx <- match(chr_df$CHROM, all_chrom)
  chr_df$Length <- as.numeric(chr_tab$Length[idx])

  seg <- data$segments
  seg <- seg[seg$CHROM %in% chr & !is.na(seg$Viterbi) &
               seg$Viterbi %in% c("D", "H"), , drop = FALSE]
  if (nrow(seg) > 0L) {
    seg <- dplyr::inner_join(
      seg,
      chr_df[, c("CHROM", "y")],
      by = "CHROM"
    )
  }

  site <- data$sites
  site <- site[site$CHROM %in% chr & !is.na(site$Genotype), , drop = FALSE]
  if (nrow(site) > 0L) {
    site <- dplyr::inner_join(
      site,
      chr_df[, c("CHROM", "y")],
      by = "CHROM"
    )
  }

  max_len <- max(chr_df$Length, na.rm = TRUE)
  x_max <- ceiling(max_len / 1e6) * 1e6
  if (x_max <= 0) {
    x_max <- max_len
  }
  breaks <- seq(0, x_max, by = 1e7)
  if (length(breaks) == 1L && x_max > 0) {
    breaks <- c(0, x_max)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = chr_df,
      mapping = ggplot2::aes(
        xmin = 0, xmax = .data$Length,
        ymin = .data$y - 1, ymax = .data$y + 1
      ),
      fill = "grey60",
      colour = "black",
      inherit.aes = FALSE
    )

  if (nrow(seg) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = seg,
      mapping = ggplot2::aes(
        xmin = .data$Start - 1, xmax = .data$End,
        ymin = .data$y, ymax = .data$y + 1,
        fill = .data$Viterbi
      ),
      colour = NA,
      inherit.aes = FALSE
    )
  }

  if (nrow(site) > 0L) {
    tick <- ggplot2::geom_segment(
      data = site,
      mapping = ggplot2::aes(
        x = .data$POS, xend = .data$POS,
        y = .data$y - 1, yend = .data$y,
        colour = .data$Genotype
      ),
      linewidth = 0.1,
      inherit.aes = FALSE
    )
    if (isTRUE(rasterize) &&
        requireNamespace("ggrastr", quietly = TRUE)) {
      tick <- ggrastr::rasterise(tick, dpi = dpi)
    }
    p <- p + tick
  }

  labels <- data$labels
  legend_labels <- c(
    D = labels$donor,
    R = labels$recurrent,
    H = labels$heterozygous
  )

  p <- p +
    ggplot2::scale_fill_manual(
      values = colors,
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = colors,
      breaks = c("R", "H", "D"),
      labels = unname(legend_labels[c("R", "H", "D")]),
      drop = FALSE,
      name = "Genotype"
    ) +
    ggplot2::scale_y_continuous(
      breaks = chr_df$y,
      labels = chr_df$Label,
      expand = ggplot2::expansion(add = 1)
    ) +
    ggplot2::scale_x_continuous(
      breaks = breaks,
      labels = paste0(breaks / 1e6, " Mb"),
      limits = c(0, x_max),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    cowplot::theme_half_open(font_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (!ext %in% c("pdf", "png")) {
      stop("`file` must have a .pdf or .png extension.", call. = FALSE)
    }
    ggplot2::ggsave(
      filename = file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  p
}

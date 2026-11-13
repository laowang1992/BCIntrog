#' Default HMM emission probability matrix
#'
#' Return the default 3 x 3 emission probability matrix used by BCIntrog.
#' Rows represent the hidden states (D, R, H) and columns represent the
#' observed genotypes (D, R, H).
#'
#' @return A 3 x 3 numeric matrix with row and column names `D`, `R` and `H`.
#'   Each row sums to 1.
#'
#' @details User-supplied matrices are not normalized. Each row must sum to 1;
#'   otherwise an error is raised.
#'
#' @examples
#' emission_probability()
#'
#' @export
emission_probability <- function() {
  matrix(
    c(
      0.90, 0.05, 0.05,
      0.05, 0.90, 0.05,
      0.40, 0.40, 0.20
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(.bt_states, .bt_states)
  )
}

.default_emission <- function() emission_probability()

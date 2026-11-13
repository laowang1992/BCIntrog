#' Physical-distance transition probability matrix
#'
#' Compute the 3 x 3 state transition probability matrix for a given physical
#' distance. The probability of at least one state change is
#' `p = 1 - exp(-distance / cross_length)`; off-diagonal transitions each have
#' probability `p / 2`.
#'
#' @param distance Physical distance between adjacent markers, in base pairs.
#'   Must be a single non-negative number.
#' @param cross_length Expected distance between state changes, in base pairs.
#'   Defaults to 10,000,000 (one change per 10 Mb).
#'
#' @return A 3 x 3 numeric matrix with row and column names `D`, `R` and `H`.
#'
#' @examples
#' transition_probability(0)
#' transition_probability(1e6, cross_length = 1e7)
#'
#' @export
transition_probability <- function(distance, cross_length = 1e7) {
  .check_scalar_numeric(distance, "distance", min = 0)
  .check_scalar_numeric(cross_length, "cross_length", min = 0)
  if (cross_length == 0) {
    stop("`cross_length` must be greater than 0.", call. = FALSE)
  }

  p <- -expm1(-distance / cross_length)
  p <- min(max(p, 0), 1 - 1e-12)

  matrix(
    c(
      1 - p, p / 2, p / 2,
      p / 2, 1 - p, p / 2,
      p / 2, p / 2, 1 - p
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(.bt_states, .bt_states)
  )
}

.transition_log_components <- function(distance, cross_length) {
  .check_scalar_numeric(cross_length, "cross_length", min = 0)
  if (cross_length == 0) {
    stop("`cross_length` must be greater than 0.", call. = FALSE)
  }
  if (anyNA(distance) || any(distance < 0)) {
    stop("Distances must be non-negative and non-missing.", call. = FALSE)
  }
  p <- -expm1(-distance / cross_length)
  p <- pmin(pmax(p, 0), 1 - 1e-12)
  list(
    log_diag = log1p(-p),
    log_off = log(p) - log(2)
  )
}

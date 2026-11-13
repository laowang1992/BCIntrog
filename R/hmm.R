#' Run a three-state HMM by Viterbi and Forward-Backward
#'
#' Run a three-state hidden Markov model (D, R, H) for one chromosome. The
#' transition probability depends on the physical distance between adjacent
#' markers. Missing observations are treated as uninformative: their emission
#' likelihood is 1 for all states.
#'
#' @param pos Numeric vector of marker positions, strictly increasing.
#' @param observations Character vector of observed genotypes (`D`, `R`, `H` or
#'   `NA`), one per marker.
#' @param emission A 3 x 3 emission probability matrix. Rows are hidden states
#'   (`D`, `R`, `H`); columns are observed states; each row must sum to 1.
#'   Defaults to [emission_probability()]. Matrices with rows not summing to 1
#'   are rejected and are not automatically normalized.
#' @param init_prob Numeric vector of initial state probabilities, length 3.
#'   Defaults to a uniform distribution.
#' @param cross_length Expected distance between state changes, in base pairs.
#'   Defaults to 10,000,000 (one change per 10 Mb).
#'
#' @return A tibble with columns `POS`, `Viterbi`, `Prob_D`, `Prob_R`, `Prob_H`
#'   and `Confidence`. If all observations are missing, `Viterbi` is `NA` and
#'   the posterior probabilities equal the initial probabilities.
#'
#' @examples
#' pos <- c(1000, 2000, 3000, 4000, 5000)
#' obs <- c("D", "D", "H", "R", "R")
#' run_hmm(pos, obs)
#'
#' @export
run_hmm <- function(pos, observations, emission = NULL, init_prob = NULL,
                    cross_length = 1e7) {
  if (!is.numeric(pos) || length(pos) == 0L) {
    stop("`pos` must be a non-empty numeric vector.", call. = FALSE)
  }
  if (anyNA(pos) || any(!is.finite(pos))) {
    stop("`pos` must not contain NA or infinite values.", call. = FALSE)
  }
  if (length(pos) > 1L && any(diff(pos) <= 0)) {
    stop("`pos` must be strictly increasing.", call. = FALSE)
  }
  if (length(observations) != length(pos)) {
    stop("`pos` and `observations` must have the same length.", call. = FALSE)
  }
  obs <- .validate_genotype(observations, "observations")
  emission <- .validate_emission(emission)
  init_prob <- .validate_init_prob(init_prob)
  .check_scalar_numeric(cross_length, "cross_length", min = 0)
  if (cross_length == 0) {
    stop("`cross_length` must be greater than 0.", call. = FALSE)
  }

  n <- length(pos)
  states <- .bt_states

  if (all(is.na(obs))) {
    warning(
      "All observations are missing; no HMM state can be inferred.",
      call. = FALSE
    )
    posterior <- matrix(
      rep(init_prob, each = n),
      nrow = n,
      ncol = 3L,
      dimnames = list(NULL, states)
    )
    return(tibble::tibble(
      POS = pos,
      Viterbi = rep(NA_character_, n),
      Prob_D = posterior[, "D"],
      Prob_R = posterior[, "R"],
      Prob_H = posterior[, "H"],
      Confidence = apply(posterior, 1L, max)
    ))
  }

  log_emission <- log(emission)
  emit <- matrix(
    0,
    nrow = n,
    ncol = 3L,
    dimnames = list(NULL, states)
  )
  for (t in seq_len(n)) {
    if (!is.na(obs[t])) {
      emit[t, ] <- log_emission[, obs[t]]
    }
  }

  if (n > 1L) {
    trans <- .transition_log_components(diff(pos), cross_length)
    log_diag <- trans$log_diag
    log_off <- trans$log_off
  }

  # Viterbi 递推：记录每个状态的最优对数概率和回溯指针
  delta <- matrix(
    -Inf,
    nrow = n,
    ncol = 3L,
    dimnames = list(NULL, states)
  )
  psi <- matrix(0L, nrow = n, ncol = 3L)
  delta[1L, ] <- log(init_prob) + emit[1L, ]

  if (n > 1L) {
    for (t in 2:n) {
      ld <- log_diag[t - 1L]
      lo <- log_off[t - 1L]
      for (j in seq_len(3L)) {
        scores <- delta[t - 1L, ] + lo
        scores[j] <- delta[t - 1L, j] + ld
        psi[t, j] <- which.max(scores)
        delta[t, j] <- max(scores) + emit[t, j]
      }
    }
  }

  if (!any(is.finite(delta[n, ]))) {
    stop(
      "No valid HMM path exists. Check the emission matrix for zero ",
      "probabilities that make all states impossible.",
      call. = FALSE
    )
  }

  path <- integer(n)
  path[n] <- which.max(delta[n, ])
  if (n > 1L) {
    for (t in (n - 1L):1L) {
      path[t] <- psi[t + 1L, path[t + 1L]]
    }
  }
  viterbi <- states[path]

  # Forward 递推（对数空间）
  alpha <- matrix(
    -Inf,
    nrow = n,
    ncol = 3L,
    dimnames = list(NULL, states)
  )
  alpha[1L, ] <- log(init_prob) + emit[1L, ]

  if (n > 1L) {
    for (t in 2:n) {
      ld <- log_diag[t - 1L]
      lo <- log_off[t - 1L]
      for (j in seq_len(3L)) {
        scores <- alpha[t - 1L, ] + lo
        scores[j] <- alpha[t - 1L, j] + ld
        alpha[t, j] <- emit[t, j] + .log_sum_exp(scores)
      }
    }
  }

  # Backward 递推（对数空间）
  beta <- matrix(
    0,
    nrow = n,
    ncol = 3L,
    dimnames = list(NULL, states)
  )
  if (n > 1L) {
    for (t in (n - 1L):1L) {
      ld <- log_diag[t]
      lo <- log_off[t]
      for (i in seq_len(3L)) {
        scores <- lo + emit[t + 1L, ] + beta[t + 1L, ]
        scores[i] <- ld + emit[t + 1L, i] + beta[t + 1L, i]
        beta[t, i] <- .log_sum_exp(scores)
      }
    }
  }

  # 后验概率按位点归一化
  posterior <- alpha + beta
  for (t in seq_len(n)) {
    denom <- .log_sum_exp(posterior[t, ])
    if (!is.finite(denom)) {
      stop(
        "Posterior probabilities could not be normalized at position ",
        pos[t],
        ". Check the emission matrix for zero probabilities.",
        call. = FALSE
      )
    }
    posterior[t, ] <- exp(posterior[t, ] - denom)
  }

  tibble::tibble(
    POS = pos,
    Viterbi = viterbi,
    Prob_D = posterior[, "D"],
    Prob_R = posterior[, "R"],
    Prob_H = posterior[, "H"],
    Confidence = apply(posterior, 1L, max)
  )
}

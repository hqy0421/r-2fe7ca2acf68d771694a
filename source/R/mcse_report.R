#' Monte Carlo error of a report computed from PSIS-LOO
#'
#' `mcse_report()` estimates the Monte Carlo standard error (MCSE) of a report
#' total that is computed by aggregating observation-wise PSIS-LOO scores. It
#' uses a chain-wise moving block bootstrap of the posterior draws to capture
#' the dependence between the pointwise scores: the same draws (and therefore
#' the same importance sampling noise) underlie every observation, so the
#' pointwise errors cannot be treated as independent and summed.
#'
#' @export
#' @param log_lik A three-dimensional array of log-likelihood values with
#'   dimensions \eqn{I} (draws per chain) by \eqn{C} (chains) by \eqn{N}
#'   (observations), as returned by [extract_log_lik()] with
#'   `merge_chains = FALSE`. All values must be finite.
#' @param weights Report aggregation weights for the \eqn{N} actual
#'   observations. Either a numeric vector of length \eqn{N} (one report
#'   total) or a numeric \eqn{N} by \eqn{R} matrix (one column per report
#'   cell, e.g. one column per region). Weights must be finite and
#'   non-negative, and each report cell must have at least one positive
#'   weight. Column names, if any, are used to name the report estimates.
#' @param block_length Length of the contiguous blocks drawn within each
#'   chain. Must be an integer between 1 and \eqn{I}. Blocks that reach the
#'   end of a chain wrap around to the beginning. If `NULL` (the default) the
#'   block length is `ceiling(sqrt(I))`, a common rule of thumb for the
#'   moving block bootstrap.
#' @param seed Integer random seed. If `NULL`, a seed is drawn from the
#'   current random number stream and returned in the result, so that the
#'   run can be reproduced. The function always uses the
#'   `"Mersenne-Twister"` / `"Inversion"` / `"Rejection"` generator and
#'   restores the caller's random number state and generator kind on exit.
#'   Reproducing the bootstrap requires the returned `seed`, `block_length`,
#'   and `n_rep`.
#' @param n_rep Number of bootstrap replicates. Must be at least 2 because
#'   the MCSE is the standard deviation of the replicate report values.
#' @param r_eff Vector of relative effective sample size estimates for the
#'   likelihood (`exp(log_lik)`) of each observation, as for [loo()]. A
#'   scalar or a vector of length \eqn{N}. The same values are used for the
#'   original data and for every bootstrap replicate, so that all report
#'   values are computed with the same PSIS tail length and remain directly
#'   comparable.
#' @template cores
#' @param k_threshold Pareto \eqn{k} threshold used to flag a replicate as
#'   unreliable. The default is the sample-size dependent threshold
#'   `1 - 1 / log10(S)` with \eqn{S = I * C}, as in
#'   [pareto-k-diagnostic]. A replicate is unreliable if any of its
#'   observation-specific \eqn{k} estimates is non-finite or exceeds the
#'   threshold. Unreliable replicates are **not** removed and are included
#'   in the MCSE.
#'
#' @details
#' For each bootstrap replicate, draw contiguous blocks of `block_length`
#' draws within each chain, with replacement, wrapping around at the end of
#' the chain, and concatenate the blocks until the original chain length is
#' reached. One resampling index per chain is applied to every observation,
#' because all observations are evaluated with the same posterior draws. The
#' resampled draws are then used to rerun PSIS and to recompute each
#' observation's LOO score. The report value for a replicate is
#' \eqn{\sum_i w_{ir} \widehat{\mathrm{elpd}}_{\mathrm{loo},i}}, using the
#' *actual* report weights (they are not resampled). The MCSE is the
#' standard deviation of these replicate report values.
#'
#' The point estimate `estimate` is computed from the original, unresampled
#' draws, so it does not depend on `seed` or `n_rep`. The bootstrap only
#' quantifies how much that estimate would move under a different posterior
#' sample.
#'
#' Replicates with high Pareto \eqn{k} are deliberately retained. Dropping
#' them would condition on the importance sampling working and would
#' understate the Monte Carlo error. The number of such replicates is
#' returned as `n_unreliable` and their positions as `unreliable`.
#'
#' The same reasoning applies to a replicate whose report value or LOO score
#' is not finite: it is counted as unreliable and kept, so the corresponding
#' `mcse` is then `NA` rather than a number computed from the remaining
#' replicates. In that case a warning with class
#' `"loo_warning_replicate_nonfinite"` is signaled.
#'
#' `mcse_report()` is separate from [mcse_loo()], which remains unchanged and
#' keeps its existing behaviour for `psis_loo` objects. In particular,
#' `mcse_loo()` combines the pointwise MCSEs analytically under the
#' independence assumption, whereas `mcse_report()` estimates the MCSE of a
#' weighted report total directly from the dependence-aware bootstrap.
#'
#' @return A named list with components:
#' \describe{
#'   \item{`estimate`}{
#'     Report total(s) computed from the original draws. A scalar, or a
#'     vector with one element per column of `weights`.
#'   }
#'   \item{`mcse`}{
#'     Bootstrap standard deviation of the replicate report values, with the
#'     same length as `estimate`.
#'   }
#'   \item{`n_unreliable`}{
#'     Number of bootstrap replicates flagged as unreliable by the Pareto
#'     \eqn{k} diagnostic or by a non-finite report value. These replicates
#'     remain in the MCSE.
#'   }
#'   \item{`unreliable`}{
#'     Logical vector of length `n_rep` flagging the unreliable replicates.
#'   }
#'   \item{`replicates`}{
#'     A `n_rep` by `R` matrix of replicate report values.
#'   }
#'   \item{`seed`, `block_length`, `n_rep`, `rng_kind`}{
#'     Parameters needed to reproduce the bootstrap.
#'   }
#'   \item{`n_chains`, `n_draws`, `n_obs`, `n_report`}{
#'     Dimensions of the input and of the report.
#'   }
#'   \item{`r_eff`, `k_threshold`}{
#'     The relative effective sample sizes and Pareto \eqn{k} threshold used.
#'   }
#'   \item{`call`}{
#'     The matched call.
#'   }
#' }
#'
#' @seealso [loo()], [mcse_loo()], [relative_eff()],
#'   [pareto-k-diagnostic]
#'
#' @examples
#' LLarr <- example_loglik_array()
#' dim(LLarr) # draws, chains, observations
#'
#' # One report total: an equally weighted mean over the 32 observations
#' w <- rep(1 / dim(LLarr)[3], dim(LLarr)[3])
#' mcse_report(LLarr, weights = w, block_length = 10, seed = 123, n_rep = 20)
#'
#' \dontrun{
#' # One report value per region (two regions here)
#' region <- rep(1:2, length.out = dim(LLarr)[3])
#' W <- cbind(region1 = as.numeric(region == 1),
#'            region2 = as.numeric(region == 2))
#' mcse_report(LLarr, weights = W, block_length = 10, seed = 123, n_rep = 200)
#' }
#'
mcse_report <- function(log_lik,
                        weights,
                        block_length = NULL,
                        seed = NULL,
                        n_rep = 1000L,
                        r_eff = 1,
                        cores = getOption("mc.cores", 1),
                        k_threshold = NULL) {
  if (!is.array(log_lik) || length(dim(log_lik)) != 3L) {
    stop(
      "'log_lik' must be a 3-D array with dimensions ",
      "(draws, chains, observations).",
      call. = FALSE
    )
  }
  if (!is.numeric(log_lik)) {
    stop("'log_lik' must be numeric.", call. = FALSE)
  }
  dims <- dim(log_lik)
  if (any(dims < 1L)) {
    stop("All dimensions of 'log_lik' must be at least 1.", call. = FALSE)
  }
  validate_log_lik(log_lik)

  I <- dims[1L] # draws per chain
  C <- dims[2L] # chains
  N <- dims[3L] # observations
  S <- I * C

  W <- .report_weights_matrix(weights, N)
  R <- ncol(W)

  if (is.null(block_length)) {
    block_length <- as.integer(ceiling(sqrt(I)))
  }
  checkmate::assert_int(block_length, lower = 1L, upper = I)
  block_length <- as.integer(block_length)

  checkmate::assert_int(n_rep, lower = 2L)
  n_rep <- as.integer(n_rep)

  if (!is.null(seed)) {
    checkmate::assert_int(seed)
    seed <- as.integer(seed)
  }

  if (is.null(k_threshold)) {
    k_threshold <- ps_khat_threshold(S)
  } else {
    checkmate::assert_number(k_threshold, finite = TRUE)
  }

  cores <- loo_cores(cores)
  r_eff <- prepare_psis_r_eff(r_eff, len = N)

  # Isolate the bootstrap from the caller's random number stream. The RNG
  # kind is pinned as well so that a stored seed reproduces the result even
  # for callers who changed RNGkind().
  rng_kind <- c("Mersenne-Twister", "Inversion", "Rejection")
  old_rng_kind <- RNGkind()
  old_has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_rng_kind)))
    if (old_has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (is.null(seed)) {
    # Draw from the caller's stream before pinning the kind, so a NULL seed
    # still depends on the caller's state in the usual way.
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  suppressWarnings(do.call(RNGkind, as.list(rng_kind)))
  set.seed(seed)

  ll_mat <- llarray_to_matrix(log_lik)

  # Point estimate from the original draws, not an average of replicates.
  # Warnings (e.g. high Pareto k or too few tail samples) are kept here so
  # that the usual loo diagnostics surface for the actual data. They are
  # muffled inside the bootstrap loop, where they would be repeated n_rep
  # times, and unreliable replicates are reported via `n_unreliable` instead.
  full <- .report_psis_elpd(ll_mat, r_eff = r_eff, cores = cores)
  estimate <- as.vector(crossprod(W, full$elpd))

  replicates <- matrix(NA_real_, nrow = n_rep, ncol = R)
  unreliable <- logical(n_rep)
  for (b in seq_len(n_rep)) {
    rows <- integer(S)
    for (cc in seq_len(C)) {
      idx <- .report_block_index(I, block_length)
      rows[((cc - 1L) * I + 1L):(cc * I)] <- (cc - 1L) * I + idx
    }
    boot <- .report_psis_elpd(
      ll_mat[rows, , drop = FALSE],
      r_eff = r_eff,
      cores = cores,
      quiet = TRUE
    )
    replicates[b, ] <- as.vector(crossprod(W, boot$elpd))
    # A replicate is unreliable if any Pareto k is non-finite or too high, or
    # if a LOO score itself did not come out finite. All replicates, reliable
    # or not, stay in the MCSE.
    bad_k <- !is.finite(boot$pareto_k) | boot$pareto_k > k_threshold
    unreliable[b] <- any(bad_k) || any(!is.finite(boot$elpd))
  }
  if (anyNA(replicates) || any(!is.finite(replicates))) {
    .warn_diagnostic(
      "Some bootstrap replicates produced non-finite report values. ",
      "The MCSE is returned as NA for the affected report cells. ",
      "See `replicates` and `n_unreliable`.",
      class = "loo_warning_replicate_nonfinite"
    )
  }

  mcse <- apply(replicates, 2L, sd)
  report_names <- colnames(W)
  if (!is.null(report_names)) {
    names(estimate) <- report_names
    names(mcse) <- report_names
    colnames(replicates) <- report_names
  }

  list(
    estimate = estimate,
    mcse = mcse,
    n_unreliable = sum(unreliable),
    unreliable = unreliable,
    replicates = replicates,
    seed = seed,
    block_length = block_length,
    n_rep = n_rep,
    rng_kind = rng_kind,
    n_chains = C,
    n_draws = I,
    n_obs = N,
    n_report = R,
    r_eff = r_eff,
    k_threshold = k_threshold,
    call = match.call()
  )
}

# internal ----------------------------------------------------------------

#' Validate and coerce report weights to an N by R matrix
#'
#' @noRd
#' @param weights Numeric vector of length N or numeric matrix with N rows.
#' @param N Number of observations.
#' @return A numeric matrix with N rows.
#'
.report_weights_matrix <- function(weights, N) {
  if (is.data.frame(weights)) {
    weights <- data.matrix(weights)
  }
  if (is.null(dim(weights))) {
    if (!is.numeric(weights)) {
      stop("'weights' must be numeric.", call. = FALSE)
    }
    if (length(weights) != N) {
      stop("'weights' must have one value per observation.", call. = FALSE)
    }
    W <- matrix(as.double(weights), ncol = 1L)
  } else {
    if (length(dim(weights)) != 2L) {
      stop("'weights' must be a vector or a matrix.", call. = FALSE)
    }
    if (!is.numeric(weights)) {
      stop("'weights' must be numeric.", call. = FALSE)
    }
    if (nrow(weights) != N) {
      stop(
        "'weights' must have one row per observation (", N, " rows).",
        call. = FALSE
      )
    }
    W <- weights
    storage.mode(W) <- "double"
  }
  if (anyNA(W) || any(!is.finite(W))) {
    stop(
      "'weights' must not contain missing or non-finite values.",
      call. = FALSE
    )
  }
  if (any(W < 0)) {
    stop("'weights' must be non-negative.", call. = FALSE)
  }
  if (any(colSums(W) <= 0)) {
    stop(
      "Each report cell in 'weights' must have at least one positive weight.",
      call. = FALSE
    )
  }
  W
}

#' Draw one chain's bootstrap indices as contiguous, wrapping blocks
#'
#' Blocks of `block_length` consecutive positions are started at uniformly
#' random positions and concatenated until `n` positions have been produced.
#' The last block is truncated. Positions wrap around the end of the chain.
#'
#' @noRd
#' @param n Number of draws in the chain.
#' @param block_length Block length.
#' @return An integer vector of length `n`.
#'
.report_block_index <- function(n, block_length) {
  n <- as.integer(n)
  block_length <- as.integer(block_length)
  n_blocks <- ceiling(n / block_length)
  starts <- sample.int(n, size = n_blocks, replace = TRUE) - 1L
  offsets <- seq.int(0L, block_length - 1L)
  idx <- rep(starts, each = block_length) + rep(offsets, times = n_blocks)
  idx <- idx %% n + 1L
  idx[seq_len(n)]
}

#' PSIS-LOO pointwise elpd and Pareto k for a log-likelihood matrix
#'
#' @noRd
#' @param ll_mat Draws by observations matrix of log-likelihood values.
#' @param r_eff Relative effective sample sizes (one per observation).
#' @param cores Number of cores for PSIS.
#' @param quiet Suppress PSIS warnings. Used for bootstrap replicates, whose
#'   reliability is summarized through the Pareto k diagnostics instead.
#' @return A list with vectors `elpd` and `pareto_k`.
#'
.report_psis_elpd <- function(ll_mat, r_eff, cores, quiet = FALSE) {
  run <- function() {
    importance_sampling.matrix(
      log_ratios = -ll_mat,
      r_eff = r_eff,
      cores = cores,
      method = "psis"
    )
  }
  psis_out <- if (quiet) suppressWarnings(run()) else run()
  lw <- weights.importance_sampling(psis_out, log = TRUE, normalize = TRUE)
  list(
    elpd = matrixStats::colLogSumExps(ll_mat + lw),
    pareto_k = psis_out$diagnostics$pareto_k
  )
}

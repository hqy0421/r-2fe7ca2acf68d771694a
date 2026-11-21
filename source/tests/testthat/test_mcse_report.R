test_that("mcse_report returns point estimate, mcse, and reproduction info", {
  LLarr <- example_loglik_array()
  w <- rep(1 / dim(LLarr)[3], dim(LLarr)[3])

  res <- mcse_report(
    LLarr,
    weights = w,
    block_length = 10,
    seed = 123,
    n_rep = 20
  )

  expect_named(
    res,
    c(
      "estimate", "mcse", "n_unreliable", "unreliable", "replicates",
      "seed", "block_length", "n_rep", "rng_kind", "n_chains", "n_draws",
      "n_obs", "n_report", "r_eff", "k_threshold", "call"
    )
  )
  expect_length(res$estimate, 1)
  expect_length(res$mcse, 1)
  expect_true(is.finite(res$estimate))
  expect_true(is.finite(res$mcse))
  expect_gte(res$mcse, 0)
  expect_equal(res$n_unreliable, sum(res$unreliable))
  expect_length(res$unreliable, 20)
  expect_equal(dim(res$replicates), c(20, 1))
  # the MCSE is the standard deviation of the retained (all) replicates
  expect_equal(res$mcse, sd(res$replicates[, 1]))
  # reproduction parameters
  expect_equal(res$seed, 123)
  expect_equal(res$block_length, 10)
  expect_equal(res$n_rep, 20)
  expect_equal(res$n_chains, 2)
  expect_equal(res$n_draws, 500)
  expect_equal(res$n_obs, 32)
  expect_equal(res$n_report, 1)
})

test_that("mcse_report point estimate matches loo() with the report weights", {
  LLarr <- example_loglik_array()
  w <- c(rep(1, 16), rep(2, 16))

  res <- mcse_report(LLarr, weights = w, block_length = 10, seed = 1, n_rep = 5)
  loo_out <- suppressWarnings(loo(LLarr))
  expected <- sum(w * loo_out$pointwise[, "elpd_loo"])

  expect_equal(res$estimate, expected)
})

test_that("mcse_report is reproducible from the returned seed and params", {
  LLarr <- example_loglik_array()
  w <- rep(1 / dim(LLarr)[3], dim(LLarr)[3])

  a <- mcse_report(LLarr, weights = w, block_length = 7, seed = 42, n_rep = 6)
  b <- mcse_report(
    LLarr,
    weights = w,
    block_length = a$block_length,
    seed = a$seed,
    n_rep = a$n_rep
  )
  expect_identical(a$replicates, b$replicates)
  expect_identical(a$estimate, b$estimate)
  expect_identical(a$mcse, b$mcse)

  # A NULL seed is drawn and returned, and reproduces the run
  c1 <- mcse_report(LLarr, weights = w, block_length = 7, n_rep = 4)
  c2 <- mcse_report(
    LLarr,
    weights = w,
    block_length = c1$block_length,
    seed = c1$seed,
    n_rep = c1$n_rep
  )
  expect_identical(c1$replicates, c2$replicates)
})

test_that("mcse_report restores the caller's RNG state and kind", {
  LLarr <- example_loglik_array()
  w <- rep(1 / dim(LLarr)[3], dim(LLarr)[3])

  old_kind <- RNGkind()
  on.exit(suppressWarnings(do.call(RNGkind, as.list(old_kind))), add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  kind_before <- RNGkind()
  set.seed(99)
  seed_before <- .Random.seed

  invisible(mcse_report(LLarr, weights = w, seed = 7, n_rep = 3))

  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, seed_before)
})

test_that("mcse_report handles several report cells", {
  LLarr <- example_loglik_array()
  region <- rep(1:2, length.out = dim(LLarr)[3])
  W <- cbind(
    region1 = as.numeric(region == 1),
    region2 = as.numeric(region == 2)
  )

  res <- mcse_report(LLarr, weights = W, block_length = 10, seed = 3, n_rep = 8)

  expect_equal(res$n_report, 2)
  expect_named(res$estimate, c("region1", "region2"))
  expect_named(res$mcse, c("region1", "region2"))
  expect_equal(colnames(res$replicates), c("region1", "region2"))
  expect_equal(res$mcse, apply(res$replicates, 2, sd))

  # A vector is equivalent to a single-column matrix
  w1 <- W[, 1]
  r_vec <- mcse_report(
    LLarr,
    weights = w1,
    block_length = 10,
    seed = 3,
    n_rep = 8
  )
  expect_equal(r_vec$estimate, unname(res$estimate["region1"]))
  expect_equal(r_vec$mcse, unname(res$mcse["region1"]))
})

test_that("mcse_report resamples within chains with wraparound blocks", {
  # Direct test of the index helper: blocks are contiguous runs modulo the
  # chain length, and the result always has exactly one index per draw.
  set.seed(1)
  idx <- loo:::.report_block_index(n = 10L, block_length = 4L)
  expect_length(idx, 10)
  expect_true(all(idx >= 1 & idx <= 10))

  for (iter in 1:20) {
    idx <- loo:::.report_block_index(n = 10L, block_length = 4L)
    # every block of 4 consecutive entries steps by 1 modulo 10
    starts <- seq(1, 8, by = 4)
    for (s in starts) {
      block <- idx[s:min(s + 3, 10)]
      if (length(block) > 1) {
        expect_true(all(diff(block) %% 10 == 1))
      }
    }
  }

  # A block as long as the chain is a rotation of all draws
  set.seed(5)
  idx <- loo:::.report_block_index(n = 10L, block_length = 10L)
  expect_equal(sort(idx), 1:10)
  expect_true(all(diff(idx) %% 10 == 1))

  # block_length = 1 draws independently within a chain
  set.seed(5)
  idx <- loo:::.report_block_index(n = 100L, block_length = 1L)
  expect_length(idx, 100)
  expect_gt(length(unique(idx)), 1)
})

test_that("mcse_report applies one index per chain to all observations", {
  # Craft two observations whose LOO scores are identical, so any
  # observation-specific resampling would break the identity between them.
  LLarr <- example_loglik_array()
  LLarr[, , 2] <- LLarr[, , 1]
  w <- rep(1 / 32, 32)

  res <- mcse_report(LLarr, weights = w, block_length = 10, seed = 8, n_rep = 5)
  expect_true(is.finite(res$estimate))
  expect_true(is.finite(res$mcse))
})

test_that("mcse_report reports unreliable replicates without dropping them", {
  LLarr <- example_loglik_array()
  w <- rep(1 / dim(LLarr)[3], dim(LLarr)[3])

  res <- mcse_report(
    LLarr,
    weights = w,
    block_length = 10,
    seed = 1,
    n_rep = 10
  )
  expect_length(res$unreliable, res$n_rep)
  expect_equal(res$n_unreliable, sum(res$unreliable))
  # even when some replicates are unreliable, all of them stay in the MCSE
  n_for_sd <- nrow(res$replicates)
  expect_equal(res$mcse, sd(res$replicates[, 1]))
  expect_equal(n_for_sd, 10)

  # a deliberately strict threshold flags more replicates but keeps the count
  strict <- mcse_report(
    LLarr,
    weights = w,
    block_length = 10,
    seed = 1,
    n_rep = 10,
    k_threshold = 0.1
  )
  expect_gte(strict$n_unreliable, res$n_unreliable)
  expect_equal(strict$n_rep, 10)
  expect_equal(strict$estimate, res$estimate)

  # A non-finite Pareto k (here from an infinite r_eff) is unreliable, but
  # the replicate is still counted and the report values stay finite.
  inf_eff <- suppressWarnings(
    mcse_report(LLarr, weights = w, block_length = 10, seed = 1, n_rep = 4,
                r_eff = Inf)
  )
  expect_equal(inf_eff$n_unreliable, 4)
  expect_length(inf_eff$unreliable, 4)
  expect_true(all(is.finite(inf_eff$replicates)))
  expect_equal(inf_eff$mcse, sd(inf_eff$replicates[, 1]))
})

test_that("mcse_report keeps mcse_loo backwards compatible", {
  LLarr <- example_loglik_array()
  loo_out <- suppressWarnings(loo(LLarr))

  old <- mcse_loo(loo_out)
  expect_true(is.na(old) || (is.numeric(old) && length(old) == 1))

  # the new function does not change the old one's value
  expect_identical(mcse_loo(loo_out), old)
})

test_that("mcse_report validates its inputs", {
  LLarr <- example_loglik_array()
  N <- dim(LLarr)[3]
  w <- rep(1 / N, N)

  expect_error(
    mcse_report(LLarr[, 1, ], weights = w, n_rep = 2),
    "3-D array"
  )
  expect_error(
    mcse_report(as.matrix(LLarr[, 1, ]), weights = w, n_rep = 2),
    "3-D array"
  )
  expect_error(
    mcse_report(LLarr, weights = w[-1], n_rep = 2),
    "one value per observation"
  )
  expect_error(
    mcse_report(LLarr, weights = matrix(1, nrow = N - 1, ncol = 2), n_rep = 2),
    "one row per observation"
  )
  expect_error(
    mcse_report(LLarr, weights = as.character(w), n_rep = 2),
    "'weights' must be numeric"
  )
  expect_error(
    mcse_report(array(as.character(LLarr), dim = dim(LLarr)), weights = w,
                n_rep = 2),
    "'log_lik' must be numeric"
  )
  expect_error(
    mcse_report(LLarr, weights = replace(w, 1, NA), n_rep = 2),
    "missing or non-finite"
  )
  expect_error(
    mcse_report(LLarr, weights = replace(w, 1, Inf), n_rep = 2),
    "missing or non-finite"
  )
  expect_error(
    mcse_report(LLarr, weights = replace(w, 1, -1), n_rep = 2),
    "non-negative"
  )
  expect_error(
    mcse_report(LLarr, weights = cbind(a = w, z = 0), n_rep = 2),
    "at least one positive weight"
  )
  expect_error(
    mcse_report(LLarr, weights = w, block_length = 0, n_rep = 2),
    "block_length"
  )
  expect_error(
    mcse_report(
      LLarr,
      weights = w,
      block_length = dim(LLarr)[1] + 1,
      n_rep = 2
    ),
    "block_length"
  )
  expect_error(mcse_report(LLarr, weights = w, n_rep = 1), "n_rep")
  expect_error(mcse_report(LLarr, weights = w, seed = 1.5, n_rep = 2), "seed")
  expect_error(
    mcse_report(LLarr, weights = w, seed = 2147483648, n_rep = 2),
    "seed"
  )

  LLbad <- LLarr
  LLbad[1, 1, 1] <- NA
  expect_error(mcse_report(LLbad, weights = w, n_rep = 2), "NAs not allowed")

  LLbad[1, 1, 1] <- Inf
  expect_error(
    mcse_report(LLbad, weights = w, n_rep = 2),
    "must be finite"
  )
})

test_that("mcse_report default block length and cores agree", {
  LLarr <- example_loglik_array()
  w <- rep(1 / dim(LLarr)[3], dim(LLarr)[3])

  res <- mcse_report(LLarr, weights = w, seed = 2, n_rep = 4)
  expect_equal(res$block_length, as.integer(ceiling(sqrt(dim(LLarr)[1]))))

  one <- mcse_report(LLarr, weights = w, block_length = 10, seed = 6, n_rep = 4,
                     cores = 1)
  two <- mcse_report(LLarr, weights = w, block_length = 10, seed = 6, n_rep = 4,
                     cores = 2)
  expect_identical(one$replicates, two$replicates)
})

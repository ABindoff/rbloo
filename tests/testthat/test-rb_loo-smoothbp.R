# Integration tests for rb_loo() on smoothbp fits. smoothbp ships a Rust
# sampler, so these are fast (no Stan compile), but they still need the package.

# The shared two-segment data generator .sbp_data() lives in helper-rb_loo.R
# (small, weakly-identified subjects -- exactly the regime where the random
# intercept drives PSIS-LOO failures); it is shared with test-reloo.R.

test_that("rb_loo cures PSIS failures on a smoothbp random-intercept fit", {
  skip_on_cran()
  skip_if_not_installed("smoothbp")

  d   <- .sbp_data()
  fit <- smoothbp::smoothbp(y ~ tau, b0 = ~ 1 + (1 | id), data = d,
                            chains = 2, iter = 1200, warmup = 600, seed = 42,
                            .verbose = FALSE)
  rb <- rb_loo(fit)

  expect_null(rb$meta$fallback)                    # RB applied, not a fallback
  expect_identical(rb$meta$model, "smoothbp")
  expect_identical(rb$meta$G, nlevels(d$id))
  expect_identical(rb$meta$N, nrow(d))

  expect_true(all(is.finite(rb$pointwise$elpd_rb)))
  expect_length(rb$structural_leverage, nrow(d))
  expect_true(all(rb$pooling_factor > 0 & rb$pooling_factor <= 1))

  # RB removes the fibre-driven failures and never adds any
  kf <- rb$diagnostics$pareto_k_full; kb <- rb$diagnostics$pareto_k
  expect_lte(sum(kb > 0.7), sum(kf > 0.7))
  if (any(kf > 0.7)) expect_lt(max(kb), max(kf))

  # structural leverage and pooling partition the group's precision exactly
  gidx <- as.integer(fit$dm$group_b0) + 1L
  for (j in unique(gidx)) {
    i <- gidx == j
    expect_equal(sum(rb$structural_leverage[i]) + rb$pooling_factor[i][1], 1,
                 tolerance = 1e-10)
  }
})

test_that("the smoothbp random-intercept downdate is exact", {
  skip_on_cran()
  skip_if_not_installed("smoothbp")

  # Conditional on a draw, the change-point mean f(t; theta) is just a number,
  # so the leave-one-out predictive marginalising u_j is the closed-form
  # Gaussian downdate. Check it against a dense numerical marginalisation.
  d   <- .sbp_data(seed = 7, J = 12)
  fit <- smoothbp::smoothbp(y ~ tau, b0 = ~ 1 + (1 | id), data = d,
                            chains = 1, iter = 600, warmup = 300, seed = 7,
                            .verbose = FALSE)
  rb  <- rb_loo(fit)

  dmat <- posterior::as_draws_matrix(fit$draws)
  sig  <- as.numeric(dmat[, "sigma"]); sigu <- as.numeric(dmat[, "sigma_u"])
  U    <- .rb_plain(dmat[, paste0("u[", fit$dm$group_levels_b0, "]"),
                         drop = FALSE])
  gidx <- as.integer(fit$dm$group_b0) + 1L
  etaF <- .rb_plain(stats::fitted(fit, summary = FALSE)) -
          U[, gidx, drop = FALSE]
  y    <- d$y

  i   <- which.max(rb$diagnostics$pareto_k_full)   # the hardest fold
  oth <- setdiff(which(gidx == gidx[i]), i)
  S   <- nrow(etaF); gr <- seq(-12, 12, length.out = 4001)
  brute <- vapply(seq_len(S), function(s) {
    uu <- sigu[s] * gr
    lw <- -0.5 * gr^2                              # log N(u | 0, sigma_u^2)
    for (o in oth) lw <- lw + dnorm(y[o], etaF[s, o] + uu, sig[s], log = TRUE)
    w <- exp(lw - max(lw))
    log(sum(w / sum(w) * dnorm(y[i], etaF[s, i] + uu, sig[s])))
  }, numeric(1))

  s2   <- sig^2
  Pmi  <- 1 / sigu^2 + length(oth) / s2
  Smi  <- rowSums(outer(rep(1, S), y[oth]) - etaF[, oth, drop = FALSE]) / s2
  clsd <- dnorm(y[i], etaF[, i] + Smi / Pmi, sqrt(s2 + 1 / Pmi), log = TRUE)
  expect_equal(clsd, brute, tolerance = 1e-8)
})

test_that("out-of-scope smoothbp fits warn and return a PSIS-LOO fallback", {
  skip_on_cran()
  skip_if_not_installed("smoothbp")

  d <- .sbp_data(seed = 3, J = 20)

  # (a) a random change-point: the mean is non-linear in the random effect
  fit_om <- smoothbp::smoothbp(y ~ tau, b0 = ~ 1 + (1 | id),
                               omega = list(~ 1 + (1 | id)), data = d,
                               chains = 1, iter = 600, warmup = 300, seed = 1,
                               .verbose = FALSE)
  expect_warning(rb_om <- rb_loo(fit_om), "random effects on omega")
  expect_psis_fallback(rb_om, "random effects on omega")

  # (b) no random intercept at all: nothing to marginalise
  fit_fe <- smoothbp::smoothbp(y ~ tau, b0 = ~ 1, data = d, chains = 1,
                               iter = 600, warmup = 300, seed = 1,
                               .verbose = FALSE)
  expect_warning(rb_fe <- rb_loo(fit_fe), "non-hierarchical")
  expect_psis_fallback(rb_fe, "non-hierarchical")

  expect_output(print(rb_om), "PSIS-LOO fallback")
})

test_that("spike-and-slab smoothbp_ss fits dispatch to the same method", {
  skip_on_cran()
  skip_if_not_installed("smoothbp")

  d   <- .sbp_data(seed = 11, J = 20)
  fit <- smoothbp::smoothbp_ss(y ~ tau, b0 = ~ 1 + (1 | id), data = d,
                               chains = 1, iter = 800, warmup = 400, seed = 11,
                               .verbose = FALSE)
  expect_s3_class(fit, "smoothbp_fit")             # inherits -> same method
  rb <- rb_loo(fit)
  expect_null(rb$meta$fallback)
  expect_true(all(is.finite(rb$pointwise$elpd_rb)))
  expect_lte(sum(rb$diagnostics$pareto_k > 0.7),
             sum(rb$diagnostics$pareto_k_full > 0.7))
})

test_that(".rb_plain strips draws_matrix so engine subsetting drops to a vector", {
  skip_if_not_installed("posterior")
  d <- posterior::as_draws_matrix(
    matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("a", "b", "c"))))
  # as.matrix() alone does NOT strip the class, and `[.draws_matrix` keeps two
  # dimensions on a single-row subset -- which breaks the engines' downdate.
  expect_s3_class(as.matrix(d[, c("a", "b"), drop = FALSE]), "draws_matrix")
  p <- .rb_plain(d[, c("a", "b"), drop = FALSE])
  expect_false(inherits(p, "draws_matrix"))
  expect_null(dim(p[3, c(1, 2)]))
  expect_length(p[3, c(1, 2)], 2L)
})

# Integration tests for rb_loo() on mcp (JAGS) fits. These need JAGS
# installed as an external program, so they are skipped wherever it is absent.

skip_no_jags <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mcp")
  testthat::skip_if_not_installed("rjags")
  testthat::skip_if(!nzchar(Sys.which("jags")) &&
                    !nzchar(Sys.which("jags.bat")),
                    "JAGS binary not on PATH")
}

# Two-segment data with a real subject-level intercept: the regime where the
# random intercept drives PSIS-LOO failures.
.mcp_data <- function(seed = 42, J = 12, nj = 8, sd_u = 2) {
  set.seed(seed)
  id  <- rep(seq_len(J), each = nj)
  x   <- rep(seq(0, 10, length.out = nj), J)
  u   <- rnorm(J, 0, sd_u)
  mu  <- 10 + u[id] + ifelse(x > 5, -1.2 * (x - 5), 0)
  data.frame(y = mu + rnorm(length(x), 0, 0.7), x = x, id = factor(id))
}

.mcp_fit <- function(d, model, seed = 1) {
  suppressWarnings(mcp::mcp(model, data = d, chains = 2, iter = 2000,
                            adapt = 1000, seed = seed))
}

test_that("rb_loo handles an mcp varying intercept exactly", {
  skip_no_jags()

  d   <- .mcp_data()
  fit <- .mcp_fit(d, list(y ~ 1 + (1 | id), 1 ~ 0 + x))
  rb  <- rb_loo(fit)

  expect_null(rb$meta$fallback)                    # RB applied, not a fallback
  expect_identical(rb$meta$model, "mcp")
  expect_identical(rb$meta$p_re, 1L)
  expect_identical(rb$meta$G, nlevels(d$id))
  expect_identical(rb$meta$N, nrow(d))
  expect_false(isTRUE(rb$meta$can_reloo))

  expect_true(all(is.finite(rb$pointwise$elpd_rb)))
  expect_length(rb$structural_leverage, nrow(d))
  expect_true(all(rb$pooling_factor > 0 & rb$pooling_factor <= 1))

  # RB removes fibre-driven failures and never adds any
  kf <- rb$diagnostics$pareto_k_full; kb <- rb$diagnostics$pareto_k
  expect_lte(sum(kb > 0.7), sum(kf > 0.7))

  # leverage and pooling partition each group's precision exactly
  gidx <- as.integer(d$id)
  for (j in unique(gidx)) {
    i <- gidx == j
    expect_equal(sum(rb$structural_leverage[i]) + rb$pooling_factor[i][1], 1,
                 tolerance = 1e-10)
  }
})

test_that("elpd_full matches mcp's own conditional loo()", {
  skip_no_jags()
  # mcp's loo() defaults to varying = TRUE, i.e. the conditional predictive.
  # rb_loo's elpd_full is PSIS over that same conditional likelihood, so the
  # two must agree; this pins the log-likelihood extraction end to end.
  #
  # This also pins the r_eff convention. .rb_engine used to pass
  # chain_id = rep(1L, S) while mcp passes the real chain ids, which changes
  # the PSIS tail length and put a 3.7e-4 gap here on a 2-chain fit. Now that
  # the engines take the real chain ids the two agree to machine precision, so
  # a regression on that convention fails this test rather than hiding in the
  # fourth decimal. See analysis/measure_reff_convention.R.
  d   <- .mcp_data(seed = 7, J = 10)
  fit <- .mcp_fit(d, list(y ~ 1 + (1 | id), 1 ~ 0 + x), seed = 7)
  rb  <- rb_loo(fit)
  lf  <- suppressWarnings(loo::loo(fit))
  expect_equal(unname(rb$estimates["elpd_full"]),
               unname(lf$estimates["elpd_loo", "Estimate"]), tolerance = 1e-8)
})

test_that("a varying change point warns and falls back to PSIS-LOO", {
  skip_no_jags()
  # The deviation shifts the break, so the mean is non-linear in u. The
  # additivity gate must catch this even though the parameter-name check
  # already would: both paths lead to the same fallback.
  d    <- .mcp_data(seed = 3, J = 10)
  cp_j <- 5 + rnorm(10, 0, 1)
  gi   <- as.integer(d$id)
  d$y  <- 10 + ifelse(d$x > cp_j[gi], -1.2 * (d$x - cp_j[gi]), 0) +
          rnorm(nrow(d), 0, 0.7)
  fit  <- .mcp_fit(d, list(y ~ 1, 1 + (1 | id) ~ 0 + x), seed = 3)

  expect_warning(rb <- rb_loo(fit), "varying change point")
  expect_psis_fallback(rb, "varying change point")
  expect_output(print(rb), "PSIS-LOO fallback")
})

test_that("a non-hierarchical mcp fit warns and falls back", {
  skip_no_jags()
  d   <- .mcp_data(seed = 5, J = 8)
  fit <- .mcp_fit(d, list(y ~ 1, 1 ~ 0 + x), seed = 5)
  expect_warning(rb <- rb_loo(fit), "non-hierarchical")
  expect_psis_fallback(rb, "non-hierarchical")
})

test_that("reloo = TRUE is refused for mcp fits, before any extraction", {
  # fit-free: the guard must fire on class alone
  expect_error(rb_loo(structure(list(), class = "mcpfit"), reloo = TRUE),
               "not implemented for mcp")
})

test_that("the additivity gate is what decides scope, not the parameter name", {
  # fit-free unit test of .rb_row_range, the gate's kernel: an additive random
  # intercept gives zero within-(draw x group) spread; anything else does not.
  delta_add <- cbind(c(1, 2), c(1, 2), c(5, 6))    # 2 draws, 3 obs
  expect_equal(max(.rb_row_range(delta_add, 1:2)), 0)
  delta_cp  <- cbind(c(1, 2), c(4, 9), c(5, 6))
  expect_gt(max(.rb_row_range(delta_cp, 1:2)), 1)
})

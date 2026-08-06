# Shared smoothbp fixture: a two-segment data set with small, weakly-identified
# subjects -- exactly the regime where the random intercept drives PSIS-LOO
# failures. Used by test-rb_loo-smoothbp.R and test-reloo.R.
.sbp_data <- function(seed = 42, J = 30) {
  set.seed(seed)
  nj  <- sample(3:5, J, replace = TRUE)
  id  <- rep(seq_len(J), nj)
  tau <- unlist(lapply(nj, function(k) sort(runif(k, 0, 10))))
  u   <- rnorm(J, 0, 2)
  mu  <- 1 + 0.4 * (tau - 5) - 1.1 * (tau - 5) / (1 + exp(-(tau - 5) * 3))
  data.frame(y = mu + u[id] + rnorm(length(id), 0, 0.5),
             tau = tau, id = factor(id))
}

# Shared expectation: a PSIS-LOO fallback object must be unmistakable -- RB
# fields all NA, PSIS fields populated, meta$fallback set. A wrong RB number is
# never returned.
expect_psis_fallback <- function(rb, reason_regex) {
  testthat::expect_s3_class(rb, "rb_loo")
  testthat::expect_false(is.null(rb$meta$fallback))
  testthat::expect_match(rb$meta$fallback, reason_regex)
  # RB-specific fields are all NA
  testthat::expect_true(all(is.na(rb$pointwise$elpd_rb)))
  testthat::expect_true(all(is.na(rb$diagnostics$pareto_k)))
  testthat::expect_true(is.na(rb$estimates[["elpd_rb"]]))
  testthat::expect_true(is.na(rb$pooling_factor))
  testthat::expect_true(is.na(rb$structural_leverage))
  # PSIS / full fields are populated and finite
  testthat::expect_true(all(is.finite(rb$pointwise$elpd_full)))
  testthat::expect_true(is.finite(rb$estimates[["elpd_full"]]))
  testthat::expect_true(is.logical(rb$refit_flag))
  testthat::expect_length(rb$diagnostics$pareto_k_full, nrow(rb$pointwise))
}

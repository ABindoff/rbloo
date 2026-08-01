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

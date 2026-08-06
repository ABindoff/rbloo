# loo_compare.rb_loo: RB-aware model comparison (all fit-free) ---------------

# a synthetic rb_loo wrapping real loo objects built from log-lik matrices
.mk_rb <- function(L, fallback = NULL) {
  lr <- suppressWarnings(loo::loo(L, r_eff = rep(1, ncol(L))))
  N  <- ncol(L)
  is_fb <- !is.null(fallback)
  structure(list(
    pooling_factor      = if (is_fb) NA_real_ else rep(0.5, N),
    structural_leverage = if (is_fb) NA_real_ else rep(0.1, N),
    pointwise = data.frame(
      elpd_rb   = if (is_fb) rep(NA_real_, N) else lr$pointwise[, "elpd_loo"],
      elpd_full = lr$pointwise[, "elpd_loo"]),
    diagnostics = list(
      pareto_k      = if (is_fb) rep(NA_real_, N) else lr$diagnostics$pareto_k,
      pareto_k_full = lr$diagnostics$pareto_k),
    refit_flag = rep(FALSE, N),
    estimates  = c(elpd_rb   = if (is_fb) NA_real_ else
                     sum(lr$pointwise[, "elpd_loo"]),
                   elpd_full = sum(lr$pointwise[, "elpd_loo"])),
    loo_rb   = if (is_fb) NA else lr,
    loo_full = lr,
    meta = list(family = "gaussian", N = N, G = 5, base_cut = 0.7,
                fallback = fallback)),
    class = "rb_loo")
}

.mk_L <- function(seed, shift = 0, S = 200, N = 40) {
  set.seed(seed)
  y <- rnorm(N)
  t(vapply(seq_len(S), function(s)
    dnorm(y, shift + rnorm(1, 0, 0.05), 1.1, log = TRUE), numeric(N)))
}

test_that("loo_compare on rb_loo objects compares the RB elpd", {
  rb1 <- .mk_rb(.mk_L(1))
  rb2 <- .mk_rb(.mk_L(2, shift = 1))   # misspecified mean: clearly worse

  cmp <- loo_compare(rb1, rb2)
  expect_s3_class(cmp, "compare.loo")
  expect_setequal(rownames(cmp), c("rb1", "rb2"))
  # the label print() actually shows (loo >= 2.9 uses a `model` column)
  if (is.data.frame(cmp)) expect_setequal(cmp$model, c("rb1", "rb2"))

  # the tabulated elpds are the RB elpds, and the ordering favours rb1
  expect_equal(cmp["rb1", "elpd_loo"], sum(rb1$pointwise$elpd_rb))
  expect_equal(cmp["rb2", "elpd_loo"], sum(rb2$pointwise$elpd_rb))
  expect_identical(rownames(cmp)[1], "rb1")
  expect_equal(cmp[1, "elpd_diff"], 0)

  # the difference and its SE match the standard paired computation
  d <- rb2$pointwise$elpd_rb - rb1$pointwise$elpd_rb
  expect_equal(unname(cmp["rb2", "elpd_diff"]), sum(d))
  expect_equal(unname(cmp["rb2", "se_diff"]), sqrt(length(d)) * sd(d))
})

test_that("mixed rb_loo / plain loo comparisons are refused", {
  rb1   <- .mk_rb(.mk_L(1))
  plain <- suppressWarnings(loo::loo(.mk_L(2), r_eff = rep(1, 40)))
  expect_error(loo_compare(rb1, plain), "different estimands")
  expect_error(loo_compare(rb1), "at least two")
})

test_that("residual flagged folds draw a warning naming the model", {
  rb1 <- .mk_rb(.mk_L(1))
  rb2 <- .mk_rb(.mk_L(2))
  rb2$refit_flag[c(3, 7)] <- TRUE
  expect_warning(loo_compare(rb1, rb2), "rb2: 2")
  expect_warning(loo_compare(rb1, rb2), "reloo = TRUE")
})

test_that("all-fallback comparisons demote to elpd_full with a warning", {
  fb1 <- .mk_rb(.mk_L(1), fallback = "unsupported family")
  fb2 <- .mk_rb(.mk_L(2, shift = 1), fallback = "unsupported family")
  expect_warning(cmp <- loo_compare(fb1, fb2), "conditional")
  expect_s3_class(cmp, "compare.loo")
  expect_equal(cmp["fb1", "elpd_loo"], sum(fb1$pointwise$elpd_full))
})

test_that("a fallback next to a real RB result is refused", {
  rb1 <- .mk_rb(.mk_L(1))
  fb2 <- .mk_rb(.mk_L(2), fallback = "unsupported family")
  expect_error(loo_compare(rb1, fb2), "fb2.*fallback")
})

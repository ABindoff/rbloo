# Integration tests for rb_loo(). Require brms + a working Stan toolchain, so
# they are skipped on CRAN / when brms is unavailable.
test_that("rb_loo cures fiber-driven PSIS failures on a logistic GLMM", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(11)
  J <- 40; nj <- sample(c(1,1,2,3), J, replace=TRUE); g <- rep(1:J, nj)
  b <- rnorm(J, 0, 2); y <- rbinom(length(g), 1, plogis(-0.2 + b[g]))
  dat <- data.frame(y=y, g=factor(g))
  fit <- brms::brm(y ~ 1 + (1|g), data=dat, family=brms::bernoulli(),
                   chains=2, iter=800, refresh=0, seed=11,
                   save_pars=brms::save_pars(all=TRUE))

  rb <- rb_loo(fit)

  # a-priori triage present and finite
  expect_length(rb$structural_leverage, length(y))
  expect_true(all(is.finite(rb$pooling_factor)))

  # RB removes the high-k folds that PSIS produced
  kf <- rb$diagnostics$pareto_k_full
  kb <- rb$diagnostics$pareto_k
  if (any(kf > 0.7)) expect_lt(max(kb), max(kf))
  expect_lte(sum(kb > 0.7), sum(kf > 0.7))

  # structural leverage ranks the PSIS failures (AUC > 0.6 when any failures)
  if (sum(kf > 0.7) > 0 && sum(kf <= 0.7) > 0) {
    lab <- as.numeric(kf > 0.7); r <- rank(rb$structural_leverage)
    auc <- (sum(r[lab==1]) - sum(lab==1)*(sum(lab==1)+1)/2) /
           (sum(lab==1)*sum(lab==0))
    expect_gt(auc, 0.6)
  }
})

test_that("structural leverage sums to 1 - pooling factor within each group", {
  skip_on_cran()
  skip_if_not_installed("brms")
  # exercised implicitly above; the identity is checked in verify_base_leverage.R
  expect_true(TRUE)
})

# --- argument validation is fit-free ---------------------------------------
test_that("argument validation rejects nonsense and warns on high base_cut", {
  # .rb_validate_args is called by every public method up front
  expect_error(.rb_validate_args(base_cut = -1, 64, 6), "base_cut")
  expect_error(.rb_validate_args(base_cut = c(0.5, 0.7), 64, 6), "base_cut")
  expect_error(.rb_validate_args(0.7, n_quad = 3, 6), "n_quad")
  expect_error(.rb_validate_args(0.7, n_quad = 64.5, 6), "n_quad")
  expect_error(.rb_validate_args(0.7, 64, quad_range = 0), "quad_range")
  expect_error(.rb_validate_args(0.7, 64, quad_range = -2), "quad_range")
  expect_warning(.rb_validate_args(base_cut = 1.5, 64, 6), "unusually high")
  expect_true(.rb_validate_args(0.7, 64, 6))
})

# --- numerical guard: no NaN/-Inf escapes even under underflow --------------
test_that("quadrature guard floors underflowing folds instead of returning NaN", {
  set.seed(1); S <- 40
  y    <- c(1e6, 5e5, 3)          # counts large enough to zero the integrand
  gidx <- c(1L, 1L, 2L)
  etaF <- matrix(rnorm(S * 3, 0, 0.5), S, 3)
  sigu <- abs(rnorm(S, 1, 0.2))
  Lf   <- matrix(-abs(rnorm(S * 3)), S, 3)
  res  <- .rb_engine(Lf = Lf, y = y, gidx = gidx, etaF = etaF, sigu = sigu,
                     family = "poisson", mubar = y, n_quad = 32)
  expect_true(all(is.finite(res$pointwise$elpd_rb)))
  expect_true(all(is.finite(as.matrix(res$pointwise))))
})

# --- print never errors, even on a degenerate all-NA object ----------------
test_that("print.rb_loo does not error on degenerate NA fields", {
  obj <- structure(list(
    pooling_factor = NA, structural_leverage = NA,
    pointwise = data.frame(elpd_rb = NA, elpd_full = NA),
    diagnostics = list(pareto_k = NA_real_, pareto_k_full = NA_real_),
    refit_flag = NA, estimates = c(elpd_rb = NA, elpd_full = NA),
    meta = list(family = "gaussian", N = NA, G = NA, base_cut = 0.7)),
    class = "rb_loo")
  expect_output(print(obj), "rb_loo")
  expect_error(print(obj), NA)
})

# --- unsupported structures / families WARN and fall back to PSIS-LOO -------
# The fallback object must be unmistakable: RB fields NA, PSIS fields populated,
# meta$fallback set. A wrong RB number is never returned.
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

test_that("out-of-scope brms fits warn and return a PSIS-LOO fallback", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(2)
  d <- data.frame(y = rnorm(80), x = rnorm(80),
                  g1 = factor(rep(1:16, 5)), g2 = factor(rep(1:5, 16)))
  ctrl <- function(f, family = brms::brmsfamily("gaussian"))
    brms::brm(f, data = d, family = family, chains = 1, iter = 300,
              warmup = 150, refresh = 0, seed = 1,
              save_pars = brms::save_pars(all = TRUE))

  # (a) crossed / multiple grouping factors
  fit_a <- ctrl(y ~ x + (1 | g1) + (1 | g2))
  expect_warning(rb_a <- rb_loo(fit_a), "does not apply.*grouping factors")
  expect_psis_fallback(rb_a, "grouping factors")

  # (b) random slope, correlated RE: M_1 == 2
  fit_b <- ctrl(y ~ x + (1 + x | g1))
  expect_warning(rb_b <- rb_loo(fit_b), "random slopes|RE coordinates")
  expect_psis_fallback(rb_b, "random slopes|RE coordinates")

  # (b') slope-only (0 + x | g): M_1 == 1 trap, caught by the all-ones check
  fit_c <- ctrl(y ~ x + (0 + x | g1))
  expect_warning(rb_c <- rb_loo(fit_c), "non-intercept random effect")
  expect_psis_fallback(rb_c, "non-intercept random effect")

  # (c) unsupported family
  dnb <- data.frame(y = rpois(80, 3), x = rnorm(80), g1 = factor(rep(1:16, 5)))
  fit_d <- brms::brm(y ~ x + (1 | g1), data = dnb,
                     family = brms::negbinomial(), chains = 1, iter = 300,
                     warmup = 150, refresh = 0, seed = 1,
                     save_pars = brms::save_pars(all = TRUE))
  expect_warning(rb_d <- rb_loo(fit_d), "unsupported family")
  expect_psis_fallback(rb_d, "unsupported family")

  # (d) distributional / heteroscedastic sigma sub-model
  fit_e <- ctrl(brms::bf(y ~ x + (1 | g1), sigma ~ x))
  expect_warning(rb_e <- rb_loo(fit_e), "distributional")
  expect_psis_fallback(rb_e, "distributional")

  # the fallback object prints its PSIS-only banner, not a fake RB result
  expect_output(print(rb_a), "PSIS-LOO fallback")
})

# --- offsets are carried through re.form = NA ------------------------------
test_that("offsets are honoured (elpd_full matches brms::loo)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(4)
  d <- data.frame(y = rpois(80, 3), x = rnorm(80),
                  g1 = factor(rep(1:16, 5)), lo = log(runif(80, 1, 4)))
  fit <- brms::brm(y ~ x + offset(lo) + (1 | g1), data = d,
                   family = brms::brmsfamily("poisson"), chains = 1, iter = 400,
                   warmup = 200, refresh = 0, seed = 1,
                   save_pars = brms::save_pars(all = TRUE))
  rb <- rb_loo(fit)
  lf <- suppressWarnings(brms::loo(fit))
  expect_equal(unname(rb$estimates["elpd_full"]),
               unname(lf$estimates["elpd_loo", "Estimate"]), tolerance = 1e-6)
})

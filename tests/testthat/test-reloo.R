# reloo = TRUE: exact refit of the residual flagged folds (brms fits only) ---

# --- argument validation and class gating are fit-free ----------------------
test_that("`reloo` is validated and rejected early on non-brms fits", {
  expect_error(.rb_validate_args(0.7, 64, 6, reloo = NA), "reloo")
  expect_error(.rb_validate_args(0.7, 64, 6, reloo = "yes"), "reloo")
  expect_error(.rb_validate_args(0.7, 64, 6, reloo = c(TRUE, FALSE)), "reloo")
  expect_true(.rb_validate_args(0.7, 64, 6, reloo = TRUE))

  # stanreg / spike-and-slab smoothbp: fail fast, BEFORE any expensive work
  expect_error(rb_loo(structure(list(), class = "stanreg"), reloo = TRUE),
               "only implemented for brmsfit")
  expect_error(
    rb_loo(structure(list(), class = c("smoothbp_ss_fit", "smoothbp_fit")),
           reloo = TRUE),
    "spike-and-slab")
})

# --- the merge is pure and unit-testable without Stan -----------------------
test_that(".rb_reloo_merge substitutes exact values and clears the flags", {
  out <- structure(list(
    pooling_factor = rep(0.5, 3), structural_leverage = rep(0.1, 3),
    pointwise   = data.frame(elpd_rb = c(-1, -5, -2), elpd_full = c(-1, -6, -2)),
    diagnostics = list(pareto_k = c(0.1, 0.9, 0.2),
                       pareto_k_full = c(0.2, 1.4, 0.3)),
    refit_flag  = c(FALSE, TRUE, FALSE),
    estimates   = c(elpd_rb = -8, elpd_full = -9),
    loo_rb = "old_loo", loo_full = "full_loo",
    meta = list(family = "gaussian", N = 3, G = 2, base_cut = 0.7)),
    class = "rb_loo")
  le <- list(
    pointwise   = matrix(c(-1, -3.5, -2), 3, 1,
                         dimnames = list(NULL, "elpd_loo")),
    diagnostics = list(pareto_k = c(0.1, 0, 0.2)))

  merged <- .rb_reloo_merge(out, le, idx = 2L)
  expect_equal(merged$pointwise$elpd_rb, c(-1, -3.5, -2))
  expect_equal(unname(merged$estimates["elpd_rb"]), -6.5)
  expect_equal(merged$diagnostics$pareto_k[2], 0)          # exact, no IS
  expect_false(any(merged$refit_flag))
  expect_identical(merged$meta$reloo_obs, 2L)
  expect_identical(merged$loo_rb, le)
  # untouched folds and the PSIS comparison are untouched
  expect_equal(merged$pointwise$elpd_full, out$pointwise$elpd_full)
  expect_equal(merged$diagnostics$pareto_k_full, out$diagnostics$pareto_k_full)
})

# --- a fallback path must not silently swallow reloo = TRUE -----------------
test_that("reloo = TRUE on a PSIS-fallback path warns that it is ignored", {
  skip_if_not_installed("brms")
  # a degenerate 'brmsfit' trips the first scope check (formula), which is
  # enough to reach fb(); the eventual fallback error is irrelevant here
  dummy <- structure(list(formula = NULL), class = "brmsfit")
  expect_warning(
    tryCatch(rb_loo(dummy, reloo = TRUE), error = function(e) NULL),
    "reloo = TRUE ignored")
})

# --- print advice is model-aware --------------------------------------------
test_that("print points brms fits at reloo = TRUE and is honest for smoothbp", {
  mk <- function(meta_extra) {
    structure(list(
      pointwise   = data.frame(elpd_rb = c(-1, -5), elpd_full = c(-1, -6)),
      diagnostics = list(pareto_k = c(0.2, 0.9), pareto_k_full = c(0.3, 1.4)),
      refit_flag  = c(FALSE, TRUE),
      estimates   = c(elpd_rb = -6, elpd_full = -7),
      meta = c(list(family = "gaussian", N = 2, G = 2, base_cut = 0.7),
               meta_extra)),
      class = "rb_loo")
  }
  expect_output(print(mk(list(can_reloo = TRUE))),
                "rb_loo\\(fit, reloo = TRUE\\)")
  expect_output(print(mk(list(model = "smoothbp", can_reloo = FALSE))),
                "update.smoothbp_fit")
  expect_output(print(mk(list())), "reloo\\(\\) or loo_moment_match\\(\\)")
  # a completed refit is reported
  done <- mk(list(can_reloo = TRUE, reloo_obs = 2L))
  done$refit_flag <- c(FALSE, FALSE)
  expect_output(print(done), "1 fold\\(s\\) refit exactly \\(obs 2\\)")
  expect_output(print(done), "no refits needed")
})

# --- end-to-end: one exact refit on a real fit (needs Stan) -----------------
test_that("rb_loo(reloo = TRUE) refits the flagged folds exactly (brms)", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(11)
  J <- 40; nj <- sample(c(1, 1, 2, 3), J, replace = TRUE); g <- rep(1:J, nj)
  b <- rnorm(J, 0, 2); y <- rbinom(length(g), 1, plogis(-0.2 + b[g]))
  dat <- data.frame(y = y, g = factor(g))
  fit <- brms::brm(y ~ 1 + (1 | g), data = dat, family = brms::bernoulli(),
                   chains = 2, iter = 800, refresh = 0, seed = 11,
                   save_pars = brms::save_pars(all = TRUE))

  rb0 <- rb_loo(fit)
  kb  <- rb0$diagnostics$pareto_k
  skip_if(max(kb) <= 0.05, "no fold with a usable RB k on this fit")

  # cut just below the worst fold: flags it (and any ties) -> minimal refits
  cut <- max(kb) - 0.01
  expect_true(isTRUE(rb0$meta$can_reloo))
  suppressWarnings(   # base_cut may exceed the 'unusually high' advisory
    rb <- rb_loo(fit, base_cut = cut, reloo = TRUE))

  i <- rb$meta$reloo_obs
  expect_identical(i, which(kb > cut))
  expect_false(any(rb$refit_flag))
  expect_true(all(rb$diagnostics$pareto_k[i] == 0))
  expect_true(all(is.finite(rb$pointwise$elpd_rb)))
  expect_true(all(rb$pointwise$elpd_rb[i] <= 0))  # bernoulli: log prob <= 0
  expect_equal(unname(rb$estimates["elpd_rb"]), sum(rb$pointwise$elpd_rb))
  # unflagged folds keep their RB values
  expect_equal(rb$pointwise$elpd_rb[-i], rb0$pointwise$elpd_rb[-i])

  # reloo = TRUE with nothing flagged is a no-op with a message, not an error
  suppressWarnings(expect_message(
    rb_hi <- rb_loo(fit, base_cut = max(kb) + 0.1, reloo = TRUE),
    "nothing to refit"))
  expect_null(rb_hi$meta$reloo_obs)
})

# --- end-to-end on smoothbp (needs a smoothbp with update()) ----------------
test_that("rb_loo(reloo = TRUE) refits flagged folds exactly (smoothbp)", {
  skip_on_cran()
  skip_if_not_installed("smoothbp")
  skip_if(is.null(utils::getS3method("update", "smoothbp_fit", optional = TRUE)),
          "installed smoothbp has no update() method")

  d   <- .sbp_data(seed = 5, J = 25)                 # shared helper fixture
  fit <- smoothbp::smoothbp(y ~ tau, b0 = ~ 1 + (1 | id), data = d,
                            chains = 2, iter = 1000, warmup = 500, seed = 5,
                            .verbose = FALSE)
  rb0 <- rb_loo(fit)
  expect_true(isTRUE(rb0$meta$can_reloo))
  kb  <- rb0$diagnostics$pareto_k
  skip_if(max(kb) <= 0.05, "no fold with a usable RB k on this fit")

  cut <- max(kb) - 0.01
  suppressWarnings(rb <- rb_loo(fit, base_cut = cut, reloo = TRUE))
  i <- rb$meta$reloo_obs
  expect_identical(i, which(kb > cut))
  expect_false(any(rb$refit_flag))
  expect_true(all(rb$diagnostics$pareto_k[i] == 0))
  expect_true(all(is.finite(rb$pointwise$elpd_rb)))
  expect_equal(rb$pointwise$elpd_rb[-i], rb0$pointwise$elpd_rb[-i])
  expect_equal(unname(rb$estimates["elpd_rb"]), sum(rb$pointwise$elpd_rb))
  # RB is exact for this scope, so refit and RB values differ only by MC noise
  expect_lt(max(abs(rb$pointwise$elpd_rb[i] - rb0$pointwise$elpd_rb[i])), 1)
})

test_that("the smoothbp singleton branch marginalises the dropped level", {
  skip_on_cran()
  skip_if_not_installed("smoothbp")
  skip_if(is.null(utils::getS3method("update", "smoothbp_fit", optional = TRUE)),
          "installed smoothbp has no update() method")

  # make group "1" a singleton so its level vanishes from the refit
  d <- .sbp_data(seed = 6, J = 12)
  d <- d[-which(d$id == "1")[-1], ]
  fit <- smoothbp::smoothbp(y ~ tau, b0 = ~ 1 + (1 | id), data = d,
                            chains = 1, iter = 800, warmup = 400, seed = 6,
                            .verbose = FALSE)
  rb   <- rb_loo(fit)
  gidx <- as.integer(fit$dm$group_b0) + 1L
  i_s  <- which(d$id == "1")
  expect_length(i_s, 1L)                             # really a singleton
  expect_identical(sum(gidx == gidx[i_s]), 1L)

  rb$refit_flag[] <- FALSE; rb$refit_flag[i_s] <- TRUE
  out <- .rb_reloo_smoothbp(rb, fit, base_cut = 0.7, gidx = gidx, y = d$y)
  expect_identical(out$meta$reloo_obs, i_s)
  expect_true(is.finite(out$pointwise$elpd_rb[i_s]))
  # for a singleton, RB-LOO's downdate is the same closed-form marginal the
  # refit branch uses, so the two agree up to MC error between posteriors
  expect_lt(abs(out$pointwise$elpd_rb[i_s] - rb$pointwise$elpd_rb[i_s]), 1)
})

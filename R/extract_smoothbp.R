# smoothbp extractor for rb_loo ---------------------------------------------
#
# smoothbp (Bindoff 2026) fits Bayesian hierarchical piecewise regression with
# logistic-smoothed change-points. The response is Gaussian and the group-level
# structure supported natively by the sampler is a random INTERCEPT on b0:
#
#     y_i | u, theta ~ N( f(t_i; theta) + u_{j(i)},  sigma^2 ),
#     u_j           ~ N( 0, sigma_u^2 ).
#
# Conditional on a draw of theta, f(t_i; theta) is just a number, so the
# leave-one-out predictive marginalising u_j is the SAME exact conjugate
# downdate used for a Gaussian random-intercept GLMM. The non-linearity of the
# change-point mean function is irrelevant to the RB algebra: it enters only
# through etaF. So RB-LOO is exact here, not an approximation.
#
# Out of scope (warn + PSIS-LOO fallback, as on the brms path):
#   * random effects on b1 / deltas / omega / rho. smoothbp codes these as
#     identity-coded dummy columns inside the design matrices (re_mask == 1),
#     shrunk by sigma_re_*. For omega/rho the mean is non-linear in u, so no
#     closed form exists; for b1/deltas the RE design column is itself a
#     function of the draw (tau - omega^(s)), which the fixed-Z engines cannot
#     represent. Either way, RB-LOO does not apply.
#   * a model with no random intercept at all (nothing to marginalise).

.rb_sbp_scope <- paste0(
  "RB-LOO scope on smoothbp fits: a random intercept on b0 -- i.e. b0 = ~ ",
  "(1 | group) -- and no random effects on b1, deltas, omega or rho.")

# Which non-b0 parameters carry random effects? smoothbp marks them with a
# 're_mask' attribute on the design matrix (1 = RE column).
.rb_sbp_re_params <- function(dm) {
  has  <- function(X) {
    m <- attr(X, "re_mask")
    !is.null(m) && any(m == 1L)
  }
  hasl <- function(L) length(L) > 0L && any(vapply(L, has, logical(1)))
  c(if (has(dm$X_b1))     "b1",
    if (hasl(dm$X_deltas)) "deltas",
    if (hasl(dm$X_om))     "omega",
    if (hasl(dm$X_rho))    "rho")
}

# Strip a posterior draws_matrix down to a plain numeric matrix. The engines
# take plain arrays: `[.draws_matrix` does NOT drop to a vector on a single-row
# subset, and the class survives arithmetic, so a draws_matrix reaching the
# engine turns etaF[s, idx] into a 1 x n_j matrix and breaks the downdate.
.rb_plain <- function(x) {
  x <- as.matrix(x)
  matrix(as.numeric(x), nrow(x), ncol(x))
}

#' @export
rb_loo.smoothbp_fit <- function(fit, base_cut = 0.7, n_quad = 64,
                                quad_range = 6, ...) {
  .rb_validate_args(base_cut, n_quad, quad_range)
  # loo::loo dispatches to smoothbp's own loo.smoothbp_fit method, so the
  # fallback needs no smoothbp-specific accessor.
  fb <- function(reason)
    .rb_psis_fallback(fit, reason, base_cut, loo_fun = function(f) loo::loo(f))

  dm <- fit$dm
  if (is.null(dm) || is.null(dm$n_groups_b0))
    stop("rb_loo: this smoothbp fit has no design-matrix block (fit$dm); it ",
         "was not produced by smoothbp()/smoothbp_ss().", call. = FALSE)

  ## ---- random effects on the non-linear / slope parameters: out of scope ----
  re_other <- .rb_sbp_re_params(dm)
  if (length(re_other) == 0L &&
      any(grepl("^sigma_re_", posterior::variables(fit$draws))))
    re_other <- "a non-b0 parameter"  # re_mask lost, but the hypers say otherwise
  if (length(re_other))
    return(fb(paste0("random effects on ", paste(re_other, collapse = " / "),
                     " (smoothbp shrinks these as design-matrix columns; for ",
                     "omega/rho the mean is non-linear in the random effect ",
                     "and for b1/deltas the effective RE design depends on the ",
                     "draw, so the RB downdate does not apply). ",
                     .rb_sbp_scope)))

  ## ---- a random intercept on b0 is required ----
  if (dm$n_groups_b0 < 1L)
    return(fb(paste0("a non-hierarchical model (no random intercept for ",
                     "RB-LOO to marginalise; plain PSIS-LOO is the usual LOO ",
                     "here). ", .rb_sbp_scope)))

  ## ---- draws: sigma, sigma_u, u[] ----
  dmat <- posterior::as_draws_matrix(fit$draws)
  vars <- colnames(dmat)
  for (v in c("sigma", "sigma_u"))
    if (!(v %in% vars))
      stop("rb_loo: smoothbp fit has no '", v, "' draw.", call. = FALSE)
  sigma <- as.numeric(dmat[, "sigma"])
  sigu  <- as.numeric(dmat[, "sigma_u"])           # an SD, not a variance
  if (any(!is.finite(sigma)) || any(sigma <= 0) ||
      any(!is.finite(sigu)) || any(sigu <= 0))
    stop("rb_loo: non-finite or non-positive 'sigma'/'sigma_u' draws.",
         call. = FALSE)

  u_nm <- paste0("u[", dm$group_levels_b0, "]")
  if (!all(u_nm %in% vars)) {
    miss <- setdiff(u_nm, vars)
    stop("rb_loo: could not match the random-intercept draws (missing: ",
         paste(miss[seq_len(min(5L, length(miss)))], collapse = ", "), ").",
         call. = FALSE)
  }
  U <- .rb_plain(dmat[, u_nm, drop = FALSE])       # S x G, group-level order

  ## ---- response and group index ----
  y <- as.numeric(fit$data[[fit$response]])
  if (is.null(y) || !length(y))
    stop("rb_loo: could not read the response '", fit$response,
         "' from fit$data.", call. = FALSE)
  if (!all(is.finite(y)))
    stop("rb_loo: response contains non-finite values.", call. = FALSE)

  gidx <- as.integer(dm$group_b0) + 1L             # smoothbp stores 0-based
  if (length(gidx) != length(y))
    stop("rb_loo: internal alignment error (response and group index differ ",
         "in length: ", length(y), " vs ", length(gidx), ").", call. = FALSE)
  if (anyNA(gidx) || any(gidx < 1L) || any(gidx > ncol(U)))
    stop("rb_loo: the grouping variable has unusable levels for some ",
         "observations (missing or out of range).", call. = FALSE)

  ## ---- fixed-effect mean (u zeroed) and the conditional log-likelihood ----
  # fitted() evaluates the full change-point mean INCLUDING u; subtracting the
  # observation's own random intercept gives the re.form = NA predictor. This
  # reuses smoothbp's own mean function rather than re-implementing it, and
  # computing Lf here avoids a second (expensive) pass through fitted().
  muF <- stats::fitted(fit, summary = FALSE)       # S x N, includes u
  if (!is.matrix(muF) || ncol(muF) != length(y) || nrow(muF) != length(sigma))
    stop("rb_loo: internal error, fitted() returned an unexpected shape (",
         paste(dim(muF), collapse = " x "), " for S=", length(sigma),
         ", N=", length(y), ").", call. = FALSE)
  muF <- .rb_plain(muF)
  Lf <- matrix(0, nrow(muF), ncol(muF))
  for (i in seq_along(y))
    Lf[, i] <- stats::dnorm(y[i], muF[, i], sigma, log = TRUE)
  etaF <- muF - U[, gidx, drop = FALSE]
  rm(muF)

  ## ---- exact Gaussian random-intercept RB-LOO (the p = 1 downdate) ----
  out <- .rb_engine(Lf = Lf, y = y, gidx = gidx, etaF = etaF, sigu = sigu,
                    family = "gaussian", sigma = sigma, trials = NULL,
                    mubar = NULL, base_cut = base_cut, n_quad = n_quad,
                    quad_range = quad_range)
  out$meta$p_re  <- 1L
  out$meta$model <- "smoothbp"
  out
}

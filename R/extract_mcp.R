# mcp (JAGS) extractor for rb_loo -------------------------------------------
#
# mcp (Lindeloev 2020) fits multiple-change-point regression via JAGS. Its
# varying effects come in two flavours, and only one of them is in RB scope:
#
#   * a varying INTERCEPT, e.g. y ~ 1 + (1 | id) on a segment. The deviation
#     enters the mean additively with a design column of ones, so conditional
#     on a draw the leave-one-out predictive marginalising u_j is the SAME
#     exact conjugate downdate used for a Gaussian random-intercept GLMM.
#     Exact, not an approximation: the change-point non-linearity enters only
#     through etaF, which we take as given.
#
#   * a varying CHANGE POINT, e.g. 1 + (1 | id) ~ 0 + x. The deviation shifts
#     the break location, so the mean is non-linear in u and no closed form
#     exists. mcp additionally draws these from a TRUNCATED normal and then
#     re-centres them (cp_1_id = cp_1_id_uncentered - mean(...)), which
#     couples the groups at O(1/G). Out of scope; PSIS-LOO fallback.
#
# Rather than infer additivity from parameter names (mcp's naming varies with
# the grouping variable and the segment index, and a varying intercept on a
# LATER segment is a jump whose design column depends on the drawn change
# point), we TEST it: fitted(varying = TRUE) - fitted(varying = FALSE) must be
# constant within each (draw, group). That admits any additive model and
# refuses any non-additive one on evidence. See analysis/probe_mcp_accessors.R,
# where the additive case measures 0 and a varying change point measures ~5.7.

.rb_mcp_scope <- paste0(
  "RB-LOO scope on mcp fits: a single varying INTERCEPT -- i.e. ",
  "y ~ 1 + (1 | group) on a segment -- with a gaussian family, no varying ",
  "change point, no varying sigma, no ARMA terms and no observation weights.")

# mcp populates fit$loglik lazily. varying = TRUE gives the CONDITIONAL
# log-likelihood (random effects included), which is the Lf the engine needs.
# NB loo.mcpfit() also defaults to varying = TRUE, so mcp's own loo() targets
# the conditional predictive -- the estimand rb_loo replaces.
.rb_mcp_loglik <- function(fit) {
  add_ll <- tryCatch(get("add_loglik", envir = asNamespace("mcp")),
    error = function(e) stop(
      "rb_loo: this version of mcp does not provide add_loglik(), which is ",
      "how the pointwise log-likelihood is populated.", call. = FALSE))
  if (is.null(fit$loglik)) fit <- add_ll(fit, varying = TRUE)
  if (is.null(fit$loglik))
    stop("rb_loo: could not obtain the pointwise log-likelihood from this mcp ",
         "fit (mcp::add_loglik() returned nothing).", call. = FALSE)
  unname(as.matrix(fit$loglik))
}

.rb_mcp_fitted <- function(fit, varying) {
  m <- stats::fitted(fit, summary = FALSE, varying = varying,
                     samples_format = "matrix")
  m <- as.matrix(m)
  matrix(as.numeric(m), nrow(m), ncol(m))
}

#' @export
rb_loo.mcpfit <- function(fit, base_cut = 0.7, n_quad = 64, quad_range = 6,
                          reloo = FALSE, ...) {
  .rb_validate_args(base_cut, n_quad, quad_range, reloo)
  if (!requireNamespace("mcp", quietly = TRUE))
    stop("rb_loo: the mcp package is needed to read this fit.", call. = FALSE)
  if (isTRUE(reloo))
    stop("rb_loo: reloo = TRUE is not implemented for mcp fits. Run ",
         "rb_loo(fit) to see which folds are flagged, then refit each ",
         "exactly by re-running mcp() without that observation.",
         call. = FALSE)

  fb <- function(reason) {
    .rb_psis_fallback(fit, reason, base_cut,
                      loo_fun = function(f) loo::loo(f))
  }

  pars <- fit$pars
  if (is.null(pars) || is.null(pars$y))
    stop("rb_loo: this mcp fit has no $pars block; it was not produced by ",
         "mcp().", call. = FALSE)

  ## ---- family / link ----
  fam  <- tryCatch(fit$family$family, error = function(e) NA_character_)
  link <- tryCatch(fit$family$link,   error = function(e) NA_character_)
  if (!identical(fam, "gaussian"))
    return(fb(paste0("family '", fam, "' on the mcp path (only gaussian is ",
                     "implemented; the non-Gaussian quadrature needs a Fisher ",
                     "weight per observation that mcp does not expose ",
                     "directly). ", .rb_mcp_scope)))
  if (!identical(link, "identity"))
    return(fb(paste0("gaussian on non-default link '", link, "'. ",
                     .rb_mcp_scope)))

  ## ---- structures the engine's conditional independence cannot survive ----
  if (length(pars$arma))
    return(fb(paste0("an ARMA/autocorrelation term (",
                     paste(pars$arma, collapse = ", "), "); the RB downdate ",
                     "assumes observations are conditionally independent ",
                     "given the random effect. ", .rb_mcp_scope)))
  if (!is.null(pars$weights))
    return(fb("observation weights (not incorporated by the RB estimator)"))
  if (length(pars$sigma) != 1L)
    return(fb(paste0("a segment-varying / distributional sigma (",
                     length(pars$sigma), " sigma parameters); the estimator ",
                     "assumes one residual scale. ", .rb_mcp_scope)))

  ## ---- exactly one varying effect, and not on cp or sigma ----
  vy <- pars$varying
  if (is.null(vy) || length(vy) == 0L)
    return(fb(paste0("a non-hierarchical model (no varying effect for RB-LOO ",
                     "to marginalise; plain PSIS-LOO is the usual LOO here). ",
                     .rb_mcp_scope)))
  if (length(vy) != 1L)
    return(fb(paste0(length(vy), " varying effects (",
                     paste(vy, collapse = ", "), "); RB-LOO handles a single ",
                     "grouping factor. ", .rb_mcp_scope)))
  if (grepl("^cp_", vy))
    return(fb(paste0("a varying change point (", vy, "): the deviation shifts ",
                     "the break location, so the mean is non-linear in the ",
                     "random effect and no closed-form downdate exists. mcp ",
                     "also truncates and re-centres these offsets, coupling ",
                     "the groups. ", .rb_mcp_scope)))
  if (grepl("^sigma", vy))
    return(fb(paste0("a varying residual scale (", vy, "). ", .rb_mcp_scope)))

  ## ---- grouping variable and index ----
  gv <- sub("^.*_[0-9]+_", "", vy)                 # Intercept_1_subj -> subj
  dat <- fit$data
  if (is.null(dat) || !(gv %in% names(dat)))
    return(fb(paste0("a varying effect whose grouping variable ('", gv,
                     "') could not be located in the fit's data")))
  gfac <- factor(dat[[gv]])
  gidx <- as.integer(gfac)
  G    <- nlevels(gfac)
  if (anyNA(gidx))
    stop("rb_loo: the grouping variable '", gv, "' has missing values.",
         call. = FALSE)

  ## ---- response ----
  y <- as.numeric(dat[[pars$y]])
  if (!length(y))
    stop("rb_loo: could not read the response '", pars$y, "' from fit$data.",
         call. = FALSE)
  if (!all(is.finite(y)))
    stop("rb_loo: response contains non-finite values.", call. = FALSE)
  N <- length(y)

  ## ---- draws: residual SD and the varying-effect SD ----
  # as_draws_*() on an mcpfit is recent; fall back to the mcmc.list, which
  # posterior:: stacks chains from in the same order fitted() does.
  dm    <- tryCatch(posterior::as_draws_matrix(fit),
                    error = function(e)
                      posterior::as_draws_matrix(fit$mcmc_post))
  dmn   <- colnames(dm)
  sd_nm <- paste0(vy, "_sd")
  if (!(pars$sigma %in% dmn))
    stop("rb_loo: mcp fit has no '", pars$sigma, "' draw.", call. = FALSE)
  if (!(sd_nm %in% dmn))
    return(fb(paste0("a varying effect whose SD parameter ('", sd_nm,
                     "') is not in the draws")))
  sigma <- as.numeric(dm[, pars$sigma])
  sigu  <- as.numeric(dm[, sd_nm])
  if (any(!is.finite(sigma)) || any(sigma <= 0) ||
      any(!is.finite(sigu))  || any(sigu <= 0))
    stop("rb_loo: non-finite or non-positive '", pars$sigma, "'/'", sd_nm,
         "' draws.", call. = FALSE)
  S <- length(sigma)

  ## ---- fixed-effect mean, and the additivity gate ----
  # fitted(varying = FALSE) is mcp's re.form = NA: the deviations are set to
  # zero (verified exactly in analysis/probe_mcp_accessors.R).
  etaF <- .rb_mcp_fitted(fit, varying = FALSE)
  f_on <- .rb_mcp_fitted(fit, varying = TRUE)
  if (!all(dim(etaF) == c(S, N)) || !all(dim(f_on) == c(S, N)))
    stop("rb_loo: internal error, mcp's fitted() returned ",
         paste(dim(f_on), collapse = " x "), " for S=", S, ", N=", N,
         ". If mcp subsampled draws, refit without an ndraws setting.",
         call. = FALSE)

  delta <- f_on - etaF                              # the RE contribution
  # Additive random intercept <=> delta is constant within (draw, group).
  tol   <- 1e-6 * max(1, stats::sd(y))
  grp   <- split(seq_len(N), gidx)
  spread <- max(vapply(grp, function(cols) {
    if (length(cols) < 2L) return(0)
    max(.rb_row_range(delta, cols))
  }, numeric(1)))
  if (!is.finite(spread) || spread > tol)
    return(fb(paste0("a varying effect that does not enter the mean ",
                     "additively (", vy, "): within a single draw and group ",
                     "the implied shift varies by up to ",
                     format(spread, digits = 3), ", where an additive random ",
                     "intercept gives 0. This is what a varying change point, ",
                     "a varying slope, or an intercept jump on a later ",
                     "segment looks like. ", .rb_mcp_scope)))

  # With additivity established, read the group effects straight off delta --
  # no need to match '<vy>[j]' column names to factor levels.
  first <- vapply(grp, function(cols) cols[1L], integer(1))
  U     <- delta[, first, drop = FALSE]             # S x G
  if (ncol(U) != G)
    stop("rb_loo: internal alignment error (", ncol(U), " group columns vs G=",
         G, ").", call. = FALSE)

  ## ---- conditional log-likelihood, with an alignment guard ----
  Lf <- .rb_mcp_loglik(fit)
  if (!all(dim(Lf) == c(S, N)))
    return(fb(paste0("a pointwise log-likelihood of shape ",
                     paste(dim(Lf), collapse = " x "), " where ", S, " x ", N,
                     " was expected")))
  # mcp is under active development with breaking changes; assert that the
  # pieces still line up rather than trusting them. dnorm(y, fitted, sigma)
  # must reproduce mcp's own loglik (exact to floating point in testing).
  chk <- unique(round(seq(1, S, length.out = min(S, 50L))))
  dev <- max(abs(vapply(seq_len(N), function(i)
    max(abs(stats::dnorm(y[i], f_on[chk, i], sigma[chk], log = TRUE) -
            Lf[chk, i])), numeric(1))))
  if (!is.finite(dev) || dev > 1e-6)
    stop("rb_loo: mcp's pointwise log-likelihood does not match ",
         "dnorm(y, fitted(varying = TRUE), ", pars$sigma, ") (max |diff| = ",
         format(dev, digits = 3), "). The draw or observation ordering of ",
         "this mcp version differs from what the extractor assumes; please ",
         "report this rather than trusting the result.", call. = FALSE)

  ## ---- exact Gaussian random-intercept RB-LOO (the p = 1 downdate) ----
  out <- .rb_engine(Lf = Lf, y = y, gidx = gidx, etaF = etaF, sigu = sigu,
                    family = "gaussian", sigma = sigma, trials = NULL,
                    mubar = NULL, base_cut = base_cut, n_quad = n_quad,
                    quad_range = quad_range)
  out$meta$p_re      <- 1L
  out$meta$model     <- "mcp"
  out$meta$can_reloo <- FALSE
  out
}

# Row-wise range over a subset of columns (per-draw spread within a group).
.rb_row_range <- function(m, cols) {
  sub <- m[, cols, drop = FALSE]
  apply(sub, 1, function(v) max(v) - min(v))
}

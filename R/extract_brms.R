# brmsfit / stanreg extractors for rb_loo -----------------------------------
#
# Design decision (maintainer): a model OUTSIDE the RB scope is NOT a hard
# error -- plain PSIS-LOO is still a valid computation for it. So every
# scope limitation (multiple/crossed/nested groups, random slopes, the
# (0+x|g) non-intercept trap, unsupported families, distributional models,
# non-default links, observation weights) emits ONE loud warning and returns
# a PSIS-LOO fallback object (see .rb_psis_fallback) that is unmistakably
# marked so it cannot be read as a real RB-LOO result. The detection logic
# below (the Z_1_1 all-ones test etc.) is unchanged; only the ACTION on a
# scope violation is warn+fallback instead of stop(). Genuine usage errors
# (a non-fit object, malformed arguments) remain hard errors.

# supported scope, for reference in the fallback reasons / docs
.rb_scope <- paste0(
  "RB-LOO scope: a SINGLE grouping factor with a random INTERCEPT only, ",
  "families gaussian / bernoulli / binomial / poisson, on their default ",
  "links (identity / logit / logit / log), with no distributional ",
  "(heteroscedastic) sub-model and no observation weights.")

#' @export
rb_loo.brmsfit <- function(fit, base_cut = 0.7, n_quad = 64, quad_range = 6, ...) {
  .rb_validate_args(base_cut, n_quad, quad_range)
  fb <- function(reason)
    .rb_psis_fallback(fit, reason, base_cut, loo_fun = brms::loo)

  ## ---- multivariate / mixture / non-linear formulas ----
  if (!inherits(fit$formula, "brmsformula"))
    return(fb("multivariate / multi-response models"))
  if (isTRUE(fit$formula$nl))
    return(fb("non-linear (nl) models"))

  ## ---- distributional / heteroscedastic sub-models ----
  # A modelled sigma (or any other dpar) breaks the single-column sigma
  # extraction below and the constant-scale assumption of the estimator.
  pforms <- fit$formula$pforms
  if (length(pforms))
    return(fb(paste0("distributional / heteroscedastic models (parameter(s) {",
                     paste(names(pforms), collapse=", "), "} are modelled)")))

  ## ---- family + link ----
  fam0 <- fit$family$family
  link <- fit$family$link
  fam  <- switch(fam0,
                 gaussian  = "gaussian",
                 bernoulli = "bernoulli",
                 binomial  = "binomial",
                 poisson   = "poisson",
                 NA_character_)
  if (is.na(fam))
    return(fb(paste0("unsupported family '", fam0, "'")))
  need_link <- switch(fam, gaussian="identity", bernoulli="logit",
                      binomial="logit", poisson="log")
  if (!identical(link, need_link))
    return(fb(paste0("family '", fam, "' on non-default link '", link,
                     "' (only '", need_link, "' is supported; the quadrature ",
                     "inverts the default link analytically)")))

  sdat <- brms::standata(fit)
  nm   <- names(sdat)

  ## ---- observation weights: the estimator ignores them ----
  if (!is.null(sdat$weights))
    return(fb("observation weights (not incorporated by the RB estimator)"))

  ## ---- exactly one grouping factor ----
  gf <- grep("^N_[0-9]+$", nm, value=TRUE)
  if (length(gf) != 1L || is.null(sdat$J_1))
    return(fb(paste0(length(gf), " grouping factors -- crossed/nested/multiple ",
                     "grouping (e.g. (1|g1)+(1|g2) or (1|g1/g2)) is out of scope")))

  ## ---- that factor must carry a single random-INTERCEPT coordinate ----
  # M_1 is the number of RE coordinates on the (only) grouping factor.
  if (is.null(sdat$M_1) || sdat$M_1 != 1L)
    return(fb(paste0("random slopes / correlated REs (", sdat$M_1,
                     " RE coordinates on the grouping factor, e.g. (1+x|g))")))
  # Z_1_1 is the RE design vector; an intercept is all-ones. A slope-only term
  # (0+x|g) also has M_1 == 1 but Z_1_1 == x, so this is the guard that stops
  # it being silently mis-handled as an intercept.
  if (is.null(sdat$Z_1_1) || !isTRUE(all(sdat$Z_1_1 == 1)))
    return(fb(paste0("a non-intercept random effect (RE design column is not ",
                     "all 1s, e.g. a random slope like (0+x|g))")))

  ## ---- smooth / spline terms (extra latent Gaussian structure) ----
  if (any(grepl("^Zs_", nm)))
    return(fb("smooth / spline terms (s(), t2(), ...)"))

  ## ---- response sanity per family (integrity checks -> hard errors) ----
  y <- as.numeric(sdat$Y)
  if (!all(is.finite(y)))
    stop("rb_loo: response contains non-finite values.", call.=FALSE)
  if (fam == "bernoulli" && !all(y %in% c(0, 1)))
    stop("rb_loo: bernoulli response must be coded 0/1; found other values.",
         call.=FALSE)
  if (fam == "poisson" && (any(y < 0) || any(y != round(y))))
    stop("rb_loo: poisson response must be non-negative integer counts.",
         call.=FALSE)
  trials <- NULL
  if (fam == "binomial") {
    if (is.null(sdat$trials))
      stop("rb_loo: binomial fit has no trials vector in standata.", call.=FALSE)
    trials <- as.numeric(sdat$trials)
    if (any(trials < 1) || any(trials != round(trials)))
      stop("rb_loo: binomial trials must be positive integers.", call.=FALSE)
    if (any(y < 0) || any(y > trials) || any(y != round(y)))
      stop("rb_loo: binomial successes must be integers in [0, trials].",
           call.=FALSE)
  }

  gidx <- as.integer(sdat$J_1)
  if (length(gidx) != length(y))
    stop("rb_loo: internal alignment error (Y and J_1 differ in length).",
         call.=FALSE)

  ## ---- draws: RE sd (exactly one) and, for gaussian, residual sd ----
  dm   <- posterior::as_draws_matrix(fit)
  vars <- posterior::variables(dm)
  sd_v <- grep("^sd_.*__Intercept$", vars, value=TRUE)
  if (length(sd_v) != 1L)
    stop("rb_loo: expected exactly one 'sd_<group>__Intercept' parameter, ",
         "found ", length(sd_v), " (", paste(sd_v, collapse=", "), "). ",
         "Cannot identify the random-intercept SD unambiguously.", call.=FALSE)
  sigu  <- as.numeric(dm[, sd_v])
  sigma <- NULL
  if (fam == "gaussian") {
    if (!("sigma" %in% vars))
      stop("rb_loo: gaussian fit has no scalar 'sigma' draw.", call.=FALSE)
    sigma <- as.numeric(dm[, "sigma"])
  }

  etaF <- brms::posterior_linpred(fit, re.form = NA)   # S x N, link scale (incl offset)
  Lf   <- brms::log_lik(fit)                           # S x N conditional
  if (!all(dim(etaF) == dim(Lf)))
    stop("rb_loo: internal error, posterior_linpred and log_lik disagree in ",
         "shape.", call.=FALSE)

  # a-priori Fisher weight per obs (posterior-mean scale)
  epred <- brms::posterior_epred(fit)                  # S x N, mean scale
  ebar  <- colMeans(epred)
  mubar <- switch(fam,
    gaussian = NULL,
    poisson  = ebar,                                   # var = mean
    bernoulli= ebar * (1 - ebar),                      # p(1-p)
    binomial = { p <- ebar / trials; trials * p * (1 - p) })

  .rb_engine(Lf=Lf, y=y, gidx=gidx, etaF=etaF, sigu=sigu, family=fam,
             sigma=sigma, trials=trials, mubar=mubar,
             base_cut=base_cut, n_quad=n_quad, quad_range=quad_range)
}

#' @export
rb_loo.stanreg <- function(fit, base_cut = 0.7, n_quad = 64, quad_range = 6, ...) {
  # rstanarm path: same estimator, rstanarm accessors.
  #
  # UNTESTED PATH. rstanarm is not installed in the package's dev environment,
  # so this method has not been exercised end-to-end. It is guarded defensively
  # and emits a one-time warning so an adopter knows they are on the
  # unvalidated path. Scope limitations warn+fall back to rstanarm::loo(),
  # exactly like the brms path.
  if (!requireNamespace("rstanarm", quietly=TRUE))
    stop("rb_loo.stanreg needs rstanarm installed.", call.=FALSE)
  .rb_validate_args(base_cut, n_quad, quad_range)
  warning("rb_loo.stanreg is an untested code path (rstanarm was not available ",
          "when the package was validated). Sanity-check the result against ",
          "brms::loo()/reloo() before relying on it.", call.=FALSE)
  fb <- function(reason)
    .rb_psis_fallback(fit, reason, base_cut, loo_fun = rstanarm::loo)

  fam0 <- fit$family$family
  link <- fit$family$link
  # rstanarm codes bernoulli as binomial with 0/1; only intercept-model those.
  fam <- switch(fam0,
                gaussian = "gaussian",
                binomial = if (max(fit$y, na.rm=TRUE) <= 1) "bernoulli"
                           else "binomial",
                poisson  = "poisson",
                NA_character_)
  if (is.na(fam))
    return(fb(paste0("unsupported family '", fam0, "'")))
  need_link <- switch(fam, gaussian="identity", bernoulli="logit",
                      binomial="logit", poisson="log")
  if (!identical(link, need_link))
    return(fb(paste0("family '", fam, "' on non-default link '", link, "'")))
  if (fam == "binomial")
    return(fb(paste0("binomial with trials > 1 on the rstanarm path (trials ",
                     "extraction is not validated)")))

  flist <- fit$glmod$reTrms$flist
  if (is.null(flist) || length(flist) != 1L)
    return(fb(paste0(length(flist), " grouping factors (crossed/nested/multiple)")))
  cnms <- fit$glmod$reTrms$cnms
  if (length(cnms) != 1L || !identical(unname(cnms[[1]]), "(Intercept)"))
    return(fb(paste0("a non-intercept random effect (term(s): ",
                     paste(unlist(cnms), collapse=", "), ")")))

  dm      <- as.matrix(fit)
  sig_col <- grep("^Sigma\\[.*Intercept.*Intercept.*\\]$", colnames(dm))
  if (length(sig_col) != 1L)
    return(fb("a random-intercept variance that could not be uniquely identified"))
  sigu  <- sqrt(dm[, sig_col])
  sigma <- NULL
  if (fam == "gaussian") {
    if (!("sigma" %in% colnames(dm)))
      stop("rb_loo.stanreg: gaussian fit has no 'sigma' column.", call.=FALSE)
    sigma <- dm[, "sigma"]
  }

  gidx <- as.integer(flist[[1]])
  etaF <- rstanarm::posterior_linpred(fit, re.form = NA)
  Lf   <- rstanarm::log_lik(fit)
  y    <- as.numeric(rstanarm::get_y(fit))
  if (fam == "bernoulli" && !all(y %in% c(0, 1)))
    stop("rb_loo.stanreg: bernoulli response must be 0/1.", call.=FALSE)
  if (fam == "poisson" && (any(y < 0) || any(y != round(y))))
    stop("rb_loo.stanreg: poisson response must be non-negative integers.",
         call.=FALSE)
  ebar  <- colMeans(rstanarm::posterior_epred(fit))
  mubar <- switch(fam, gaussian=NULL, poisson=ebar, bernoulli=ebar*(1-ebar))

  .rb_engine(Lf=Lf, y=y, gidx=gidx, etaF=etaF, sigu=sigu, family=fam,
             sigma=sigma, trials=NULL, mubar=mubar, base_cut=base_cut,
             n_quad=n_quad, quad_range=quad_range)
}

# =====================================================================
# rbloo :: Rao-Blackwellised leave-one-out for hierarchical models
#
# rb_loo(fit) returns, per observation, an a-priori pooling-factor / structural
# leverage triage (Gelman-Pardoe / fibr), the RB-LOO pointwise elpd with its
# base Pareto-k diagnostic, and a refit flag for the residual folds. RB-LOO
# marginalises the random-effect block and importance-samples only the base, so
# the fiber-driven PSIS-LOO failures disappear at zero refit cost. See Bindoff
# (2026, fibr) for the pooling factor / structural leverage, and Vehtari, Gelman
# & Gabry (2017) for PSIS-LOO and the k-hat diagnostic.
#
# Scope, by family, on a single grouping factor:
#   * gaussian  -- ANY random-effect design (intercept and/or slopes, correlated
#                  or not). The conditional is exactly Gaussian, so RB-LOO is the
#                  closed-form p x p matrix downdate (.rb_engine_gaussian_mv):
#                  EXACT, no approximation.
#   * bernoulli / binomial / poisson -- a single random INTERCEPT, via cheap 1-D
#                  quadrature over the RE (.rb_engine): numerically exact to grid
#                  resolution. Random slopes on these families need multivariate
#                  quadrature (not yet implemented) and warn + fall back to PSIS.
# =====================================================================

#' Rao-Blackwellised LOO
#'
#' @param fit a fitted model (brmsfit; stanreg via the shared extractor)
#' @param base_cut Pareto-k threshold above which a fold is flagged for refit
#' @param n_quad number of quadrature nodes for non-Gaussian families
#' @param quad_range half-width of the standardised RE quadrature grid, in prior
#'   standard deviations. The default (6) is ample for grouped REs; for
#'   observation-level (singleton) REs with extreme responses the integrand can
#'   sit further out, so widen to 8-10.
#' @return an object of class "rb_loo"
#' @export
rb_loo <- function(fit, base_cut = 0.7, n_quad = 64, quad_range = 6, ...)
  UseMethod("rb_loo")

#' @export
rb_loo.default <- function(fit, ...)
  stop("rb_loo: no method for class ", paste(class(fit), collapse="/"),
       ". Supported: brmsfit, stanreg.", call.=FALSE)

# ---------------------------------------------------------------------
# argument validation (shared by all methods). Fails loudly and early,
# before any expensive extraction, with actionable messages.
# ---------------------------------------------------------------------
.rb_validate_args <- function(base_cut, n_quad, quad_range) {
  if (!is.numeric(base_cut) || length(base_cut) != 1L ||
      !is.finite(base_cut) || base_cut <= 0)
    stop("rb_loo: `base_cut` must be a single positive finite number ",
         "(a Pareto-k refit threshold, typically 0.5-0.7).", call.=FALSE)
  if (base_cut >= 1)
    warning("rb_loo: `base_cut` = ", base_cut, " is unusually high for a ",
            "Pareto-k threshold (typically 0.5-0.7); few/no folds will be flagged.",
            call.=FALSE)
  if (!is.numeric(n_quad) || length(n_quad) != 1L || !is.finite(n_quad) ||
      n_quad != round(n_quad) || n_quad < 8)
    stop("rb_loo: `n_quad` must be a single integer >= 8 ",
         "(number of quadrature nodes; default 64).", call.=FALSE)
  if (!is.numeric(quad_range) || length(quad_range) != 1L ||
      !is.finite(quad_range) || quad_range <= 0)
    stop("rb_loo: `quad_range` must be a single positive finite number ",
         "(RE grid half-width in prior SDs; default 6).", call.=FALSE)
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# PSIS-LOO fallback for models OUTSIDE the RB scope. RB-LOO cannot be
# computed, but plain PSIS-LOO is a valid (if potentially unreliable)
# computation for these fits, so instead of erroring we warn loudly and
# return an rb_loo object that is UNMISTAKABLY marked as PSIS-only: the
# RB-specific fields (pooling_factor, structural_leverage, elpd_rb,
# pareto_k) are all NA, only the PSIS/full fields are populated, and
# meta$fallback carries the reason. print.rb_loo surfaces this prominently.
# ---------------------------------------------------------------------
.rb_advice <- paste0(
  "PSIS-LOO can be unreliable on influential observations -- inspect the ",
  "Pareto-k diagnostics, and for folds with k>0.7 use reloo() or ",
  "loo_moment_match() to correct them.")

.rb_psis_fallback <- function(fit, reason, base_cut, loo_fun) {
  warning("RB-LOO does not apply to ", reason, "; returning plain PSIS-LOO. ",
          .rb_advice, call.=FALSE)
  lf <- suppressWarnings(loo_fun(fit))
  ek <- lf$pointwise[, "elpd_loo"]
  kk <- lf$diagnostics$pareto_k
  N  <- length(ek)
  fam <- tryCatch(fit$family$family, error=function(e) NA_character_)
  if (is.null(fam) || length(fam) != 1L) fam <- NA_character_
  structure(list(
    pooling_factor      = NA_real_,
    structural_leverage = NA_real_,
    pointwise           = data.frame(elpd_rb = rep(NA_real_, N), elpd_full = ek),
    diagnostics         = list(pareto_k = rep(NA_real_, N), pareto_k_full = kk),
    refit_flag          = kk > base_cut,
    estimates           = c(elpd_rb = NA_real_, elpd_full = sum(ek)),
    loo_rb = NA, loo_full = lf,
    meta = list(family = fam, N = N, G = NA, base_cut = base_cut,
                fallback = reason)
  ), class = "rb_loo")
}

# ---------------------------------------------------------------------
# core engine. All model-specific extraction is done by the S3 methods,
# which hand the engine plain arrays:
#   Lf    : S x N full conditional log-lik matrix   (from log_lik(fit))
#   y     : response (N)
#   gidx  : group index 1..G per obs (N)
#   etaF  : S x N fixed-effect linear predictor, link scale (re.form = NA)
#   sigu  : S draws of the RE standard deviation
#   sigma : S draws of residual SD           (gaussian only)
#   trials: N binomial trials                (binomial only)
#   mubar : N a-priori Fisher weight per obs (see .fisher_weight)
# ---------------------------------------------------------------------
.rb_engine <- function(Lf, y, gidx, etaF, sigu, family, sigma=NULL, trials=NULL,
                       mubar=NULL, base_cut=0.7, n_quad=64, quad_range=6) {
  # args are validated up front by the public S3 methods (.rb_validate_args)
  S <- nrow(etaF); N <- ncol(etaF); G <- max(gidx)
  loo_of <- function(L) {
    reff <- tryCatch(loo::relative_eff(exp(L), chain_id=rep(1L,S)),
                     error=function(e) rep(1,N))
    suppressWarnings(loo::loo(L, r_eff=reff))
  }

  ## ---- RB-LOO: marginalise the RE block given each base draw ----
  Lrb <- matrix(0, S, N)
  grid <- seq(-quad_range, quad_range, length.out = n_quad)   # standardised RE nodes
  lp0  <- -0.5 * grid^2                             # log N(u|0,sigu^2) up to const & Jacobian
  for (i in 1:N) {
    j <- gidx[i]; oth <- which(gidx==j); oth <- oth[oth!=i]
    if (family == "gaussian") {
      s2  <- sigma^2
      Smi <- if (length(oth))
                rowSums(outer(rep(1,S), y[oth]) - etaF[,oth,drop=FALSE]) / s2 else rep(0,S)
      Pmi <- 1/sigu^2 + length(oth)/s2
      Lrb[,i] <- dnorm(y[i], etaF[,i] + Smi/Pmi, sqrt(s2 + 1/Pmi), log=TRUE)
    } else {
      node <- outer(sigu, grid)                    # S x Q  : u = sigu * standardised node
      lp   <- matrix(lp0, S, n_quad, byrow=TRUE)
      for (o in oth) lp <- lp + .llk(y[o], etaF[,o] + node, family, trials[o])
      # self-normalised importance weights over the RE grid, guarded against
      # underflow (all-node -Inf rows -> uniform fallback) so no NaN escapes.
      mx  <- apply(lp, 1, max)
      fin <- is.finite(mx)
      w   <- exp(lp - ifelse(fin, mx, 0))
      rs  <- rowSums(w)
      ok  <- fin & is.finite(rs) & rs > 0
      if (any(ok))  w[ok, ]  <- w[ok, ] / rs[ok]
      if (any(!ok)) w[!ok, ] <- 1 / n_quad         # degenerate draw: fall back to prior
      pred <- rowSums(w * exp(.llk(y[i], etaF[,i] + node, family, trials[i])))
      Lrb[,i] <- log(pmax(pred, 1e-300))
    }
  }

  lf <- loo_of(Lf); lr <- loo_of(Lrb)
  k_full <- lf$diagnostics$pareto_k; e_full <- lf$pointwise[,"elpd_loo"]
  k_base <- lr$diagnostics$pareto_k; e_rb   <- lr$pointwise[,"elpd_loo"]

  ## ---- a-priori pooling factor + structural leverage (fibr) ----
  su_h <- mean(sigu)
  w_i  <- if (family=="gaussian") rep(1/mean(sigma)^2, N) else mubar
  Wj   <- vapply(1:G, function(j) sum(w_i[gidx==j]) + 1/su_h^2, numeric(1))
  h_struct <- w_i / Wj[gidx]
  pi_j     <- (1/su_h^2) / Wj

  structure(list(
    pooling_factor      = pi_j[gidx],
    structural_leverage = h_struct,
    pointwise           = data.frame(elpd_rb=e_rb, elpd_full=e_full),
    diagnostics         = list(pareto_k=k_base, pareto_k_full=k_full),
    refit_flag          = (k_base > base_cut),
    estimates           = c(elpd_rb=sum(e_rb), elpd_full=sum(e_full)),
    loo_rb=lr, loo_full=lf,
    meta = list(family=family, N=N, G=G, base_cut=base_cut)
  ), class="rb_loo")
}

# =====================================================================
# Exact multivariate Gaussian RB-LOO.
#
# For a group with a p-vector of random effects u_j ~ N(0, Sigma) and RE design
# Z_j (n_j x p), the full conditional given hypers is exactly Gaussian:
#     Lambda_j = Sigma^{-1} + Z_j' Z_j / sigma^2      (p x p posterior precision)
#     V_j      = Lambda_j^{-1},   m_j = V_j Z_j' r_j / sigma^2
# Leaving out observation i downdates by rank one (Sherman-Morrison):
#     Lambda_j^{-i} = Lambda_j - z_i z_i' / sigma^2,  b_j^{-i} = b_j - z_i r_i/sigma^2
# and the leave-i-out predictive is exactly
#     y_i ~ N( eta_i + z_i' m_j^{-i},  sigma^2 + z_i' V_j^{-i} z_i ).
# No quadrature, no linearisation: this is exact for correlated random effects
# (random intercept + slopes). The scalar random-intercept downdate in
# .rb_engine is the p = 1 special case.
#
# Arrays: Z (N x p) RE design; Sig (S x p x p) draws of the RE covariance;
# sigma (S) residual SD; etaF (S x N) fixed predictor; Lf (S x N) conditional
# log-lik for the PSIS comparison. Cost is O(S * G) small p x p solves plus
# O(S * N) rank-1 updates.
# =====================================================================
.rb_engine_gaussian_mv <- function(Lf, y, gidx, etaF, Z, Sig, sigma, base_cut=0.7) {
  S <- nrow(etaF); N <- ncol(etaF); p <- ncol(Z); G <- max(gidx)
  loo_of <- function(L) {
    reff <- tryCatch(loo::relative_eff(exp(L), chain_id=rep(1L,S)),
                     error=function(e) rep(1,N))
    suppressWarnings(loo::loo(L, r_eff=reff))
  }
  grp_idx <- split(seq_len(N), gidx)

  ## ---- RB-LOO: exact matrix downdate per group and draw ----
  Lrb <- matrix(0, S, N)
  for (j in seq_len(G)) {
    idx <- grp_idx[[j]]; if (is.null(idx)) next
    Zj  <- Z[idx, , drop=FALSE]; nj <- length(idx); ZtZ <- crossprod(Zj)
    for (s in seq_len(S)) {
      s2   <- sigma[s]^2
      Vinv <- solve(Sig[s, , ]) + ZtZ / s2          # Lambda_j
      V    <- solve(Vinv)                            # V_j (p x p)
      r    <- y[idx] - etaF[s, idx]                  # residuals (n_j)
      b    <- crossprod(Zj, r) / s2                  # Z_j' r_j / s2  (p)
      for (ii in seq_len(nj)) {
        zi    <- Zj[ii, ]
        Vz    <- as.numeric(V %*% zi)
        denom <- 1 - sum(zi * Vz) / s2               # > 0 always (see notes)
        Vmi   <- V + tcrossprod(Vz) / (s2 * denom)   # V_j^{-i} (Sherman-Morrison)
        mmi   <- as.numeric(Vmi %*% (b - zi * (r[ii] / s2)))
        pm    <- etaF[s, idx[ii]] + sum(zi * mmi)
        pv    <- s2 + sum(zi * (Vmi %*% zi))
        Lrb[s, idx[ii]] <- dnorm(y[idx[ii]], pm, sqrt(pv), log=TRUE)
      }
    }
  }

  lf <- loo_of(Lf); lr <- loo_of(Lrb)
  k_full <- lf$diagnostics$pareto_k; e_full <- lf$pointwise[,"elpd_loo"]
  k_base <- lr$diagnostics$pareto_k; e_rb   <- lr$pointwise[,"elpd_loo"]

  ## ---- a-priori generalized leverage + pooling (posterior mean) ----
  # h_i = z_i' V_j z_i / sigma^2  (generalized hat value; reduces to the scalar
  # structural leverage when p = 1). pi_j = tr(Sigma^{-1} V_j)/p in (0,1],
  # the multivariate Gelman-Pardoe pooling factor (1 = fully pooled).
  Sbar <- apply(Sig, c(2, 3), mean); s2b <- mean(sigma)^2
  SbarInv <- solve(Sbar)
  h_struct <- numeric(N); pi_grp <- numeric(G)
  for (j in seq_len(G)) {
    idx <- grp_idx[[j]]; if (is.null(idx)) next
    Zj  <- Z[idx, , drop=FALSE]
    Vb  <- solve(SbarInv + crossprod(Zj) / s2b)
    for (ii in seq_along(idx)) h_struct[idx[ii]] <- sum(Zj[ii,] * (Vb %*% Zj[ii,])) / s2b
    pi_grp[j] <- sum(diag(SbarInv %*% Vb)) / p
  }

  structure(list(
    pooling_factor      = pi_grp[gidx],
    structural_leverage = h_struct,
    pointwise           = data.frame(elpd_rb=e_rb, elpd_full=e_full),
    diagnostics         = list(pareto_k=k_base, pareto_k_full=k_full),
    refit_flag          = (k_base > base_cut),
    estimates           = c(elpd_rb=sum(e_rb), elpd_full=sum(e_full)),
    loo_rb=lr, loo_full=lf,
    meta = list(family="gaussian", N=N, G=G, base_cut=base_cut, p_re=p)
  ), class="rb_loo")
}

# per-family conditional log-likelihood on a matrix of linear predictors
.llk <- function(yv, eta, family, trials=NULL)
  switch(family,
    bernoulli = if (yv==1) plogis(eta, log.p=TRUE) else plogis(-eta, log.p=TRUE),
    binomial  = dbinom(yv, trials, plogis(eta), log=TRUE),
    poisson   = dpois(yv, exp(eta), log=TRUE),
    stop("unsupported family in .llk: ", family))

#' @export
print.rb_loo <- function(x, ...) {
  m <- x$meta
  # defensive helpers: never error on NA / empty / degenerate fields
  n_gt   <- function(v, t) sum(v > t, na.rm=TRUE)
  safemax <- function(v) { v <- v[is.finite(v)]; if (length(v)) max(v) else NA_real_ }
  fnum   <- function(v) if (length(v) && is.finite(v)) sprintf("%.1f", v) else "NA"
  fk     <- function(v) if (is.finite(v)) sprintf("%.2f", v) else "NA"
  fG     <- function(v) if (is.null(v) || length(v) != 1 || is.na(v)) "NA" else as.character(v)
  kf <- x$diagnostics$pareto_k_full; kb <- x$diagnostics$pareto_k

  ## ---- PSIS-LOO fallback: RB-LOO was NOT applied ----
  if (!is.null(m$fallback)) {
    cat("rb_loo  [PSIS-LOO fallback -- RB-LOO NOT applied]\n")
    cat(sprintf("  reason  : %s\n", m$fallback))
    cat(sprintf("  elpd (PSIS-LOO) = %s    (N=%s obs)\n",
                fnum(x$estimates["elpd_full"]), fG(m$N)))
    cat(sprintf("  PSIS-LOO : #(k>0.7) = %d  (max %s)\n",
                n_gt(kf, 0.7), fk(safemax(kf))))
    cat("  advice  : PSIS-LOO can be unreliable on influential observations;\n")
    cat("            for folds with k>0.7 use reloo() or loo_moment_match().\n")
    return(invisible(x))
  }

  cat(sprintf("rb_loo  (%s GLMM; N=%s obs, G=%s groups)\n",
              m$family, fG(m$N), fG(m$G)))
  cat(sprintf("  elpd_rb = %s    elpd_full(PSIS) = %s\n",
              fnum(x$estimates["elpd_rb"]), fnum(x$estimates["elpd_full"])))
  cat(sprintf("  PSIS-LOO : #(k>0.7) = %d  (max %s)\n", n_gt(kf, 0.7), fk(safemax(kf))))
  cat(sprintf("  RB-LOO   : #(k>0.7) = %d  (max %s)   <- fiber failures removed\n",
              n_gt(kb, 0.7), fk(safemax(kb))))
  cat(sprintf("  refit_flag: %d fold(s) still k_base>%.2f -> send to reloo\n",
              sum(x$refit_flag, na.rm=TRUE), m$base_cut))
  invisible(x)
}

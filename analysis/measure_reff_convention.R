# =====================================================================
# How much does .rb_engine's r_eff convention actually change?
#
# .rb_engine calls loo::relative_eff(exp(L), chain_id = rep(1L, S)), i.e. it
# treats all draws as ONE chain. Every downstream package passes the real
# chain ids. This script measures the consequence on the paper's real-data
# case (E-real: brms::epilepsy, Poisson, patient random intercept, 4 chains)
# for BOTH log-likelihood matrices, and on every quantity the paper reports.
#
# Self-validating: it rebuilds L_rb independently and checks that the
# reconstruction reproduces rb_loo()'s own elpd_rb before comparing anything.
#
# Structural leverage and the pooling factor are computed from Fisher blocks
# and posterior means, NOT from importance weights, so they are invariant to
# r_eff by construction. Only k-hat and elpd can move.
#
# Prints summaries only. No rows.
# =====================================================================
suppressMessages({
  library(brms); library(loo); library(posterior)
  devtools::load_all("C:/Users/bindoffa/antigravity_projects/fibr_sbc",
                     quiet = TRUE)
})
options(mc.cores = 4)

hdr <- function(x) cat("\n== ", x, " ", strrep("=", max(0, 56 - nchar(x))), "\n", sep = "")

hdr("fit (E-real, as analysis/ereal_epilepsy.R)")
data(epilepsy, package = "brms")
epilepsy$zBase <- as.numeric(scale(log(epilepsy$Base)))
epilepsy$zAge  <- as.numeric(scale(log(epilepsy$Age)))
fit <- brm(count ~ zBase * Trt + zAge + (1 | patient), data = epilepsy,
           family = poisson(), chains = 4, iter = 1500, refresh = 0, seed = 7,
           save_pars = save_pars(all = TRUE))

rb   <- rb_loo(fit)
y    <- epilepsy$count
gidx <- as.integer(factor(epilepsy$patient))
N    <- length(y); G <- max(gidx)

L_full <- brms::log_lik(fit)
S      <- nrow(L_full)
nch    <- brms::nchains(fit)
cid    <- rep(seq_len(nch), each = S / nch)          # brms stacks by chain
cat(sprintf("  N = %d, G = %d, S = %d, chains = %d\n", N, G, S, nch))

## ---- rebuild L_rb independently (1-D quadrature, poisson intercept) ----
dm    <- as_draws_matrix(fit)
sdv   <- grep("^sd_.*__Intercept$", variables(fit), value = TRUE)
sigu  <- as.numeric(dm[, sdv])
etaF  <- brms::posterior_linpred(fit, re.form = NA)
gb    <- seq(-6, 6, length.out = 64)
L_rb  <- matrix(0, S, N)
for (i in seq_len(N)) {
  j <- gidx[i]; oth <- setdiff(which(gidx == j), i)
  node <- outer(sigu, gb)
  lp   <- matrix(-0.5 * gb^2, S, length(gb), byrow = TRUE)
  for (o in oth) lp <- lp + dpois(y[o], exp(etaF[, o] + node), log = TRUE)
  mx <- apply(lp, 1, max); w <- exp(lp - mx); w <- w / rowSums(w)
  L_rb[, i] <- log(pmax(rowSums(w * dpois(y[i], exp(etaF[, i] + node))), 1e-300))
}

## ---- the three conventions in the codebase ----
psis <- function(L, mode) {
  reff <- switch(mode,
    one   = loo::relative_eff(exp(L), chain_id = rep(1L, nrow(L))),  # .rb_engine
    real  = loo::relative_eff(exp(L), chain_id = cid),               # downstream
    unity = rep(1, ncol(L)))                                         # some scripts
  suppressWarnings(loo::loo(L, r_eff = reff))
}
modes <- c("one", "real", "unity")
Lf <- lapply(modes, function(m) psis(L_full, m)); names(Lf) <- modes
Lr <- lapply(modes, function(m) psis(L_rb,   m)); names(Lr) <- modes

hdr("self-check: reconstruction matches rb_loo()")
cat(sprintf("  elpd_full  rebuilt vs rb_loo : %.6f vs %.6f  (diff %.2e)\n",
            unname(Lf$one$estimates["elpd_loo", "Estimate"]),
            unname(rb$estimates["elpd_full"]),
            unname(Lf$one$estimates["elpd_loo", "Estimate"] -
                   rb$estimates["elpd_full"])))
cat(sprintf("  elpd_rb    rebuilt vs rb_loo : %.4f vs %.4f  (diff %.2e)\n",
            unname(Lr$one$estimates["elpd_loo", "Estimate"]),
            unname(rb$estimates["elpd_rb"]),
            unname(Lr$one$estimates["elpd_loo", "Estimate"] -
                   rb$estimates["elpd_rb"])))
cat("  (elpd_rb is quadrature-sensitive; a small gap here is the grid, not r_eff)\n")

## ---- what actually moves ----
cmp <- function(A, B, lab) {
  ea <- unname(A$estimates["elpd_loo", "Estimate"])
  eb <- unname(B$estimates["elpd_loo", "Estimate"])
  ka <- A$diagnostics$pareto_k; kb <- B$diagnostics$pareto_k
  cross <- sum((ka > 0.7) != (kb > 0.7))
  cat(sprintf("  %-28s d(elpd) = %+8.5f | max|dk| = %.5f | folds crossing 0.7: %d\n",
              lab, eb - ea, max(abs(ka - kb)), cross))
  invisible(cross)
}
hdr("conditional log-lik (drives k_full: the AUC/Spearman target)")
c1 <- cmp(Lf$one, Lf$real,  "chain_id=1 -> real chains")
c2 <- cmp(Lf$one, Lf$unity, "chain_id=1 -> r_eff=1")
hdr("RB log-lik (drives elpd_rb and the 'cured' claim)")
cmp(Lr$one, Lr$real,  "chain_id=1 -> real chains")
cmp(Lr$one, Lr$unity, "chain_id=1 -> r_eff=1")

## ---- the paper's reported statistics, recomputed under each convention ----
h   <- rb$structural_leverage                      # r_eff-invariant
sp  <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
auc <- function(s, l) {
  if (sum(l) == 0 || sum(l) == length(l)) return(NA_real_)
  r <- rank(s); (sum(r[l == 1]) - sum(l == 1) * (sum(l == 1) + 1) / 2) /
    (sum(l == 1) * sum(l == 0))
}
hdr("paper statistics under each convention")
cat("  convention   #k_full>0.7   Spearman(h,k_full)   AUC(h ranks k>0.7)\n")
for (m in modes) {
  k <- Lf[[m]]$diagnostics$pareto_k
  cat(sprintf("  %-11s  %11d   %18.4f   %17s\n", m, sum(k > 0.7), sp(h, k),
              formatC(auc(h, as.numeric(k > 0.7)), format = "f", digits = 4)))
}
hdr("cured claim (# RB folds still k>0.7) under each convention")
for (m in modes)
  cat(sprintf("  %-11s  %d\n", m, sum(Lr[[m]]$diagnostics$pareto_k > 0.7)))

hdr("verdict")
cat(sprintf("  Threshold crossings on k_full: %d (real chains), %d (r_eff=1).\n",
            c1, c2))
cat("  Zero crossings => the reported failure counts, Spearman and AUC are\n")
cat("  unchanged, and the fix is cosmetic for this case. Nonzero => the\n")
cat("  reported numbers move and the experiments need re-running.\n")

# =====================================================================
# Validate the multivariate-GLMM RB-LOO (tier 2): a Poisson random intercept +
# slope. RB-LOO (p-D quadrature) should track brute-force refit fold-by-fold,
# far better than PSIS-LOO, in the small-group failing regime.
# =====================================================================
suppressMessages({library(brms); library(loo); library(matrixStats)})
options(mc.cores=4)
for (f in list.files("../R", full.names=TRUE)) source(f)   # run from analysis/

set.seed(21)
J <- 30; nj <- sample(c(2,2,3), J, replace=TRUE); N <- sum(nj); g <- rep(1:J, nj)
x <- rnorm(N)
a <- rnorm(J, 0, 1.0); b <- rnorm(J, 0, 0.8)               # random intercept + slope
eta <- -0.3 + a[g] + (0.4 + b[g]) * x
y <- rpois(N, exp(eta))
dat <- data.frame(y=y, x=x, g=factor(g))
cat(sprintf("Poisson random-slope GLMM: J=%d groups, N=%d obs\n", J, N))

fit <- brm(y ~ x + (1 + x | g), data=dat, family=poisson(),
           chains=4, iter=1500, refresh=0, seed=21, backend="rstan",
           save_pars=save_pars(all=TRUE), control=list(adapt_delta=0.95))

rb <- rb_loo(fit)
stopifnot(is.null(rb$meta$fallback), rb$meta$p_re == 2)     # MV engine ran
cat(sprintf("rb_loo MV-GLMM engine: p_re=%d\n", rb$meta$p_re))
k_full <- rb$diagnostics$pareto_k_full; k_base <- rb$diagnostics$pareto_k
e_full <- rb$pointwise$elpd_full;       e_rb   <- rb$pointwise$elpd_rb
cat(sprintf("PSIS-LOO #(k>0.7)=%d (max %.2f);  RB-LOO #(k>0.7)=%d (max %.2f)\n",
            sum(k_full>0.7), max(k_full), sum(k_base>0.7), max(k_base)))

## ---- brute-force exact refit on a stratified sample of folds ----
set.seed(3)
samp <- unique(c(order(-k_full)[1:8], sample(which(k_full<0.5), 8)))
exact <- sapply(samp, function(i){
  fi <- update(fit, newdata=dat[-i,], refresh=0, recompile=FALSE)
  ll <- log_lik(fi, newdata=dat[i,,drop=FALSE])
  logSumExp(ll) - log(length(ll))
})
D <- data.frame(k=k_full[samp], e_rb=e_rb[samp], e_full=e_full[samp], e_exact=exact)
rmse <- function(a,b) sqrt(mean((a-b)^2))
cat("\n================ MV-GLMM RB-LOO exactness ================\n")
cat(sprintf("  folds refit: %d   (k range %.2f-%.2f)\n", nrow(D), min(D$k), max(D$k)))
cat(sprintf("  RMSE vs exact refit:  RB-LOO=%.4f   PSIS-LOO=%.4f\n",
            rmse(D$e_rb,D$e_exact), rmse(D$e_full,D$e_exact)))
cat(sprintf("  mean|error|:          RB-LOO=%.4f   PSIS-LOO=%.4f\n",
            mean(abs(D$e_rb-D$e_exact)), mean(abs(D$e_full-D$e_exact))))
hi <- D$k>0.7
if (any(hi)) cat(sprintf("  high-k folds only:    RB-LOO=%.4f   PSIS-LOO=%.4f\n",
            rmse(D$e_rb[hi],D$e_exact[hi]), rmse(D$e_full[hi],D$e_exact[hi])))
cat("\n  per-fold (sorted by k):\n")
print(round(D[order(D$k), c("k","e_rb","e_full","e_exact")], 4), row.names=FALSE)

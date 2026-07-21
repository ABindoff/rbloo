# =====================================================================
# Validate the exact multivariate-Gaussian RB-LOO: a random intercept+slope
# LMM. RB-LOO should match brute-force refit fold-by-fold (exact), while
# PSIS-LOO drifts on the influential folds.
# =====================================================================
suppressMessages({library(brms); library(loo); library(matrixStats)})
options(mc.cores=4)
for (f in list.files("../R", full.names=TRUE)) source(f)   # run from analysis/

set.seed(11)
J <- 25; nj <- sample(c(2,2,3), J, replace=TRUE); N <- sum(nj); g <- rep(1:J, nj)
x <- rnorm(N)
Sig <- matrix(c(0.9^2, 0.4*0.9*0.9, 0.4*0.9*0.9, 0.9^2), 2, 2)
ab <- t(chol(Sig)) %*% matrix(rnorm(2*J), 2, J)            # (a_j, b_j)
y <- 0.5 + ab[1,g] + (1.0 + ab[2,g])*x + rnorm(N, 0, 0.7)
dat <- data.frame(y=y, x=x, g=factor(g))
cat(sprintf("random-slope LMM: J=%d groups, N=%d obs\n", J, N))

fit <- brm(y ~ x + (1 + x | g), data=dat, family=gaussian(),
           chains=4, iter=1500, refresh=0, seed=11, backend="rstan",
           save_pars=save_pars(all=TRUE), control=list(adapt_delta=0.95))

rb <- rb_loo(fit)
stopifnot(is.null(rb$meta$fallback))                        # must NOT be a fallback
cat(sprintf("rb_loo ran the MV engine: p_re=%d, fallback=%s\n",
            rb$meta$p_re, !is.null(rb$meta$fallback)))
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
cat("\n================ MV-Gaussian RB-LOO exactness ================\n")
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

## ---- regression check: plain (1|g) gaussian still works (p_re=1) ----
dat2 <- dat; fit2 <- update(fit, formula. = y ~ x + (1|g), newdata=dat2, refresh=0)
rb2 <- rb_loo(fit2)
cat(sprintf("\n(1|g) gaussian: p_re=%d, fallback=%s, RB #(k>0.7)=%d (was PSIS %d)\n",
            rb2$meta$p_re, is.null(rb2$meta$fallback),
            sum(rb2$diagnostics$pareto_k>0.7), sum(rb2$diagnostics$pareto_k_full>0.7)))

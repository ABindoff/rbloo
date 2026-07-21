# Decisive bounded test: on the OLRE model, does the EXACT refit track RB (-651.6)
# or PSIS/reloo-0.7 (-643.7)? Refit a stratified sample of folds across the k range
# -- crucially including the 0.4<k<0.7 folds that reloo(0.7) leaves at their PSIS
# value -- and compare RB vs PSIS vs exact refit per fold.
suppressMessages({library(brms); library(loo); library(posterior)})
options(mc.cores=4)
pkg_R <- if (dir.exists("../R")) "../R" else "R"   # package lives at repo root
for (f in list.files(pkg_R, pattern="\\.R$", full.names=TRUE)) source(f)

data(epilepsy, package="brms")
epilepsy$zBase <- as.numeric(scale(log(epilepsy$Base)))
epilepsy$zAge  <- as.numeric(scale(log(epilepsy$Age)))
epilepsy$obs   <- factor(seq_len(nrow(epilepsy)))
N <- nrow(epilepsy)

M1 <- brm(count ~ zBase*Trt + zAge + (1|obs), data=epilepsy, family=poisson(),
          chains=4, iter=2000, refresh=0, seed=7, backend="rstan",
          save_pars=save_pars(all=TRUE), control=list(adapt_delta=0.99, max_treedepth=12))
lw <- suppressWarnings(loo(M1)); k <- lw$diagnostics$pareto_k
e_psis <- lw$pointwise[,"elpd_loo"]
rb <- rb_loo(M1, n_quad=96, quad_range=8); e_rb <- rb$pointwise$elpd_rb

set.seed(3)
# stratified sample across k buckets, emphasising the disagreement zone 0.4-0.7
buck <- cut(k, c(-Inf,0.4,0.55,0.7,Inf))
samp <- unlist(lapply(split(seq_len(N), buck), function(ix) sample(ix, min(length(ix), 10))))
cat(sprintf("refitting %d folds across k buckets: %s\n", length(samp),
            paste(names(table(buck[samp])), table(buck[samp]), sep=":", collapse="  ")))

# exact leave-one-out refit for each sampled fold (marginal predictive)
exact <- sapply(samp, function(i){
  fi <- update(M1, newdata=epilepsy[-i,], refresh=0, recompile=FALSE)
  ll <- log_lik(fi, newdata=epilepsy[i,,drop=FALSE], allow_new_levels=TRUE,
                sample_new_levels="gaussian")
  matrixStats::logSumExp(ll) - log(length(ll))
})

D <- data.frame(obs=samp, k=k[samp], e_psis=e_psis[samp], e_rb=e_rb[samp], e_exact=exact)
D$err_psis <- D$e_psis - D$e_exact; D$err_rb <- D$e_rb - D$e_exact
rmse <- function(x) sqrt(mean(x^2))
cat("\n================ ADJUDICATION: exact refit vs RB vs PSIS ================\n")
for(b in levels(buck)){ ii<-cut(D$k,c(-Inf,0.4,0.55,0.7,Inf))==b; if(!any(ii)) next
  cat(sprintf("  k in %-12s (n=%2d):  mean(PSIS-exact)=%+.3f  mean(RB-exact)=%+.3f\n",
              b, sum(ii), mean(D$err_psis[ii]), mean(D$err_rb[ii]))) }
cat(sprintf("\n  ALL sampled folds:  RMSE(PSIS-exact)=%.3f   RMSE(RB-exact)=%.3f\n", rmse(D$err_psis), rmse(D$err_rb)))
cat(sprintf("  mean signed bias :  PSIS %+.3f (over-optimistic if >0)   RB %+.3f\n", mean(D$err_psis), mean(D$err_rb)))
cat(sprintf("  sub-0.7 folds only: RMSE(PSIS)=%.3f  RMSE(RB)=%.3f   <- the folds reloo(0.7) does NOT refit\n",
            rmse(D$err_psis[D$k<0.7]), rmse(D$err_rb[D$k<0.7])))
saveRDS(D, "adjudicate.rds"); print(round(D[order(D$k),c("k","e_psis","e_rb","e_exact","err_psis","err_rb")],3), row.names=FALSE)

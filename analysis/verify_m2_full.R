# Airtight M2: correct the negative-binomial model's LOO to exact (reloo at a low
# threshold, since NB refits are cheap), so BOTH models in E-decision use an
# exact-quality elpd. Then recompute the verdict ladder consistently.
suppressMessages({library(brms); library(loo)})
options(mc.cores=4)

data(epilepsy, package="brms")
epilepsy$zBase <- as.numeric(scale(log(epilepsy$Base)))
epilepsy$zAge  <- as.numeric(scale(log(epilepsy$Age)))
N <- nrow(epilepsy)

M2 <- brm(count ~ zBase*Trt + zAge, data=epilepsy, family=negbinomial(),
          chains=4, iter=2000, refresh=0, seed=7, backend="rstan",
          save_pars=save_pars(all=TRUE), control=list(adapt_delta=0.95))
lw <- suppressWarnings(loo(M2)); k <- lw$diagnostics$pareto_k
e_psis2 <- lw$pointwise[,"elpd_loo"]
cat(sprintf("M2 (NB): PSIS total %.1f  max k=%.2f  #(k>0.2)=%d\n", sum(e_psis2), max(k), sum(k>0.2)))

# correct EVERY fold where PSIS could be biased (k>0.2); NB refits are fast
t0<-Sys.time()
rl <- suppressWarnings(reloo(M2, loo=lw, k_threshold=0.2))
e_corr2 <- rl$pointwise[,"elpd_loo"]
cat(sprintf("M2 corrected (reloo k>0.2, %d refits, %.0fs): total %.1f  (moved %+.1f nats)\n",
            sum(k>0.2), as.numeric(Sys.time()-t0,units="secs"), sum(e_corr2), sum(e_corr2)-sum(e_psis2)))

# M1 vectors (validated exact via RB, RMSE 0.017): reuse saved OLRE results
z1 <- readRDS("ereal_epilepsy_olre.rds")
e_psis1 <- z1$e_full; e_rb1 <- z1$e_rb; e_gold1 <- z1$e_gold

cmp <- function(e1,e2,lab){ d<-e1-e2; De<-sum(d); SE<-sqrt(length(d))*sd(d); z<-De/SE
  cat(sprintf("  %-34s  elpd_diff = %+6.1f  (SE %4.1f, z=%+.1f)  [%s]\n", lab, De, SE, z,
              ifelse(abs(z)>=2,"significant","not significant"))); c(De=De,SE=SE) }

cat("\n================ E-decision, with M2 corrected to exact ================\n")
cat(sprintf("M1 total elpd:  PSIS %.1f | RB %.1f | reloo(0.7) %.1f\n", sum(e_psis1),sum(e_rb1),sum(e_gold1)))
cat(sprintf("M2 total elpd:  PSIS %.1f | corrected(reloo k>0.2) %.1f\n\n", sum(e_psis2), sum(e_corr2)))
r_psis <- cmp(e_psis1, e_psis2, "PSIS-LOO (naive, both models)")
r_relo <- cmp(e_gold1, e_corr2, "reloo k>0.7 (M1) vs exact M2")
r_rb   <- cmp(e_rb1,   e_corr2, "RB-LOO (M1) vs exact M2  <-- honest")
saveRDS(list(r_psis=r_psis,r_relo=r_relo,r_rb=r_rb,
             e_psis1=e_psis1,e_rb1=e_rb1,e_gold1=e_gold1,e_psis2=e_psis2,e_corr2=e_corr2),
        "e_decision.rds")
cat(sprintf("\nM2 correction = %+.1f nats (%s vs the model gap)\n",
            sum(e_corr2)-sum(e_psis2),
            ifelse(abs(sum(e_corr2)-sum(e_psis2))<2,"negligible","modest")))

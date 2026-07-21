# =====================================================================
# E-decision : does PSIS-LOO's over-optimism flip a real model-selection choice?
#
# Two standard ways to model overdispersed counts, compared on brms::epilepsy:
#   M1  count ~ zBase*Trt + zAge + (1|obs)   Poisson-lognormal (observation-level RE)
#   M2  count ~ zBase*Trt + zAge             negative binomial (no latent per obs)
# These are near-equivalent honest predictive models, but PSIS-LOO fails on 97 of
# M1's folds (the per-obs RE) and is over-optimistic there, while M2's PSIS is
# well-behaved -> PSIS spuriously prefers M1. RB-LOO (cheap) reproduces the exact
# refit (reloo) decision; PSIS does not.
#
# M1's expensive vectors (PSIS / RB / reloo pointwise elpd) are reused from
# ereal_epilepsy_olre.rds (the 82-min reloo already run). We fit M2 here.
# =====================================================================
suppressMessages({library(brms); library(loo); library(posterior)})
options(mc.cores=4)

data(epilepsy, package="brms")
epilepsy$zBase <- as.numeric(scale(log(epilepsy$Base)))
epilepsy$zAge  <- as.numeric(scale(log(epilepsy$Age)))
epilepsy$obs   <- factor(seq_len(nrow(epilepsy)))
N <- nrow(epilepsy)

## ---- M1 (OLRE Poisson): reuse saved pointwise elpd vectors (data order) ----
z1 <- readRDS("ereal_epilepsy_olre.rds")
stopifnot(length(z1$e_full)==N, all(z1$count==epilepsy$count))   # alignment check
e_psis1 <- z1$e_full     # PSIS-LOO (biased on the 97 failing folds)
e_rb1   <- z1$e_rb       # RB-LOO   (marginalises the per-obs RE, cheap)
e_gold1 <- z1$e_gold     # reloo    (exact refit of the 97 folds; 82 min)

## ---- M2 (negative binomial, no RE): fit + loo + reloo ----
M2 <- brm(count ~ zBase*Trt + zAge, data=epilepsy, family=negbinomial(),
          chains=4, iter=2000, refresh=0, seed=7, backend="rstan",
          save_pars=save_pars(all=TRUE), control=list(adapt_delta=0.95))
loo2 <- suppressWarnings(loo(M2))
cat(sprintf("M2 (NB): PSIS #(k>0.7)=%d / %d  max k=%.2f\n",
            sum(loo2$diagnostics$pareto_k>0.7), N, max(loo2$diagnostics$pareto_k)))
reloo2 <- suppressWarnings(reloo(M2, loo=loo2, k_threshold=0.7))
e_psis2 <- loo2$pointwise[,"elpd_loo"]
e_gold2 <- reloo2$pointwise[,"elpd_loo"]   # honest M2 (few refits, cheap)

## ---- model comparison: elpd difference (M1 - M2) and its SE ----
# SE of the pointwise-elpd difference = sqrt(N) * sd(diff)  (standard loo_compare SE)
compare <- function(e1, e2, label){
  d <- e1 - e2; Delta <- sum(d); SE <- sqrt(length(d))*sd(d)
  z <- Delta/SE
  pick <- if(Delta>0) "M1 (OLRE Poisson)" else "M2 (neg-binomial)"
  sig  <- if(abs(z)>=2) "SIGNIFICANT" else "not significant"
  cat(sprintf("  %-28s  elpd_diff(M1-M2) = %+7.1f  (SE %4.1f, z=%+.1f)  -> prefer %-18s [%s]\n",
              label, Delta, SE, z, pick, sig))
  c(Delta=Delta, SE=SE)
}

cat("\n================ E-decision: model selection, OLRE-Poisson vs neg-binomial ================\n")
cat(sprintf("M1 total elpd:  PSIS %.1f | RB %.1f | reloo %.1f\n", sum(e_psis1), sum(e_rb1), sum(e_gold1)))
cat(sprintf("M2 total elpd:  PSIS %.1f | reloo %.1f\n\n", sum(e_psis2), sum(e_gold2)))
cat("Model-selection verdict under each estimator:\n")
r_psis <- compare(e_psis1, e_psis2, "PSIS-LOO (naive)")          # what a practitioner gets by default
r_rb   <- compare(e_rb1,   e_gold2, "RB-LOO (+reloo on NB)")     # cheap recommended
r_gold <- compare(e_gold1, e_gold2, "reloo exact (gold)")        # truth

flip <- sign(r_psis["Delta"]) != sign(r_gold["Delta"]) ||
        (abs(r_psis["Delta"]/r_psis["SE"])>=2) != (abs(r_gold["Delta"]/r_gold["SE"])>=2)
cat(sprintf("\n>>> PSIS decision matches gold? %s   |  RB decision matches gold? %s\n",
            ifelse(!flip,"yes","NO -- PSIS FLIPS THE DECISION"),
            ifelse(sign(r_rb["Delta"])==sign(r_gold["Delta"]) &&
                   (abs(r_rb["Delta"]/r_rb["SE"])>=2)==(abs(r_gold["Delta"]/r_gold["SE"])>=2),
                   "yes","no")))

saveRDS(list(e_psis1=e_psis1,e_rb1=e_rb1,e_gold1=e_gold1,e_psis2=e_psis2,e_gold2=e_gold2,
             r_psis=r_psis,r_rb=r_rb,r_gold=r_gold), "e_decision.rds")

## ---- figure: forest plot of elpd_diff +/- 2 SE under the three estimators ----
png("e_decision.png", width=1100, height=520, res=135)
par(mar=c(4.6,9,3,2))
labs <- c("reloo exact\n(gold; 82 min)","RB-LOO\n(cheap)","PSIS-LOO\n(naive)")
D  <- c(r_gold["Delta"], r_rb["Delta"], r_psis["Delta"])
SE <- c(r_gold["SE"],    r_rb["SE"],    r_psis["SE"])
cols <- c("#2ca02c","#1f77b4","#d62728")
xlim <- range(c(D-2*SE, D+2*SE, 0)) + c(-3,3)
plot(D, 1:3, pch=19, col=cols, cex=1.6, yaxt="n", ylab="", xlim=xlim, ylim=c(0.5,3.5),
     xlab="elpd difference  (M1 OLRE-Poisson  -  M2 neg-binomial)",
     main="Same data, same two models:\nPSIS-LOO picks M1; the exact answer says tie")
segments(D-2*SE, 1:3, D+2*SE, 1:3, col=cols, lwd=3)
abline(v=0, lty=2, col="grey50")
axis(2, at=1:3, labels=labs, las=1, cex.axis=0.85)
text(D, 1:3+0.22, sprintf("%+.1f +/- %.1f", D, 2*SE), cex=0.8, col=cols)
mtext("M2 better  <---            --->  M1 better", side=1, line=3, cex=0.8, col="grey40")
dev.off(); cat("\nWrote figure: e_decision.png\n")

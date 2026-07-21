# =====================================================================
# E-real (OLRE stressor): epilepsy with an observation-level random effect.
#   count ~ zBase*Trt + zAge + (1|obs)   -- Poisson-lognormal overdispersion.
# Every observation is its own singleton group -> the RE is maximally
# data-driven -> PSIS-LOO fails hard on the extreme counts. RB-LOO marginalises
# the per-obs RE (1-D quadrature over the lognormal) = the integrated predictive.
# Single grouping factor, so rb_loo() applies directly (quad_range widened for
# the singleton tails). Head-to-head vs reloo-exact gold + moment matching.
# =====================================================================
suppressMessages({library(brms); library(loo); library(posterior)})
options(mc.cores=4)
pkg_R <- if (dir.exists("../R")) "../R" else "R"   # package lives at repo root
for (f in list.files(pkg_R, pattern="\\.R$", full.names=TRUE)) source(f)

data(epilepsy, package="brms")
epilepsy$zBase <- as.numeric(scale(log(epilepsy$Base)))
epilepsy$zAge  <- as.numeric(scale(log(epilepsy$Age)))
epilepsy$obs   <- factor(seq_len(nrow(epilepsy)))

fit <- brm(count ~ zBase*Trt + zAge + (1|obs), data=epilepsy,
           family=poisson(), chains=4, iter=2000, refresh=0, seed=7,
           backend="rstan", save_pars=save_pars(all=TRUE),
           control=list(adapt_delta=0.99, max_treedepth=12))

## ---- rb_loo() : widened quadrature range for the singleton tails ----
rb <- rb_loo(fit, n_quad=96, quad_range=8)
print(rb)

k_full <- rb$diagnostics$pareto_k_full; k_base <- rb$diagnostics$pareto_k
e_full <- rb$pointwise$elpd_full;       e_rb   <- rb$pointwise$elpd_rb
h_struct <- rb$structural_leverage
hi <- which(k_full>0.7); N <- length(k_full)
cat(sprintf("\nPSIS-LOO: #(k>0.7)=%d / %d  max k=%.2f\n", length(hi), N, max(k_full)))

sp<-function(a,b) suppressWarnings(cor(a,b,method="spearman"))
auc<-function(s,l){ if(sum(l)==0||sum(l)==length(l)) return(NA); r<-rank(s)
  (sum(r[l==1])-sum(l==1)*(sum(l==1)+1)/2)/(sum(l==1)*sum(l==0)) }
rmse<-function(a,b) sqrt(mean((a-b)^2))

lw <- suppressWarnings(loo(fit))
t0<-Sys.time(); rl <- suppressWarnings(reloo(fit, loo=lw, k_threshold=0.7))
reloo_secs<-as.numeric(Sys.time()-t0,units="secs"); e_gold<-rl$pointwise[,"elpd_loo"]

t0<-Sys.time(); mm <- tryCatch(suppressWarnings(loo_moment_match(fit, loo=lw, k_threshold=0.7)),
                               error=function(e){cat("MM ERR:",conditionMessage(e),"\n");NULL})
mm_secs<-as.numeric(Sys.time()-t0,units="secs")
e_mm <- if(!is.null(mm)) mm$pointwise[,"elpd_loo"] else rep(NA,N)
k_mm <- if(!is.null(mm)) mm$diagnostics$pareto_k else rep(NA,N)

cat("\n============ E-real OLRE : epilepsy Poisson-lognormal overdispersion ============\n")
cat(sprintf("C1  Spearman(structural leverage, k_full)=%.2f   AUC(k>0.7)=%.2f\n",
            sp(h_struct,k_full), auc(h_struct,as.numeric(k_full>0.7))))
cat(sprintf("C2  cure: PSIS #(k>0.7)=%d -> RB #(k>0.7)=%d  (mean k %.2f -> %.2f; max %.2f -> %.2f)\n",
            sum(k_full>0.7), sum(k_base>0.7), mean(k_full), mean(k_base), max(k_full), max(k_base)))
cat(sprintf("    MM #(k>0.7 after)=%d\n", sum(k_mm>0.7, na.rm=TRUE)))
cat(sprintf("Acc vs reloo-exact on %d high-k folds (elpd RMSE): PSIS=%.3f  MM=%.3f  RB=%.3f\n",
            length(hi), rmse(e_full[hi],e_gold[hi]), rmse(e_mm[hi],e_gold[hi]), rmse(e_rb[hi],e_gold[hi])))
cat(sprintf("Cost: RB refits/opts=0 | MM opts=%d (%.0fs) | reloo refits=%d (%.0fs)\n",
            length(hi), mm_secs, length(hi), reloo_secs))

saveRDS(list(k_full=k_full,k_base=k_base,k_mm=k_mm,e_full=e_full,e_rb=e_rb,e_mm=e_mm,
             e_gold=e_gold,h_struct=h_struct,hi=hi,pooling=rb$pooling_factor,
             count=epilepsy$count),
        "ereal_epilepsy_olre.rds")

## ---- figure ----
png("ereal_epilepsy_olre.png", width=1600, height=520, res=135)
par(mfrow=c(1,3), mar=c(4.3,4.4,3.2,1))
plot(h_struct, k_full, pch=19, col=adjustcolor("#333333",0.5),
     xlab="structural leverage (a-priori)", ylab="PSIS-LOO k-hat",
     main=sprintf("C1: predicts k-hat\nrho=%.2f AUC=%.2f",
                  sp(h_struct,k_full), auc(h_struct,as.numeric(k_full>0.7))))
abline(h=0.7,lty=3,col="red")
plot(k_full, k_base, pch=19, col="#2ca02c",
     xlab="k-hat full (PSIS)", ylab="k-hat base (RB-LOO)",
     main=sprintf("C2: RB cures %d failures", length(hi)))
abline(0,1,lty=2,col="grey60"); abline(h=0.7,v=0.7,lty=3,col="grey70")
yl<-range(c(e_gold[hi],e_rb[hi],e_full[hi],e_mm[hi]),na.rm=TRUE)
plot(e_gold[hi], e_rb[hi], pch=19, col="#1f77b4", xlim=yl, ylim=yl,
     xlab="reloo exact elpd (gold)", ylab="LOO elpd",
     main=sprintf("accuracy (%d high-k folds)",length(hi)))
points(e_gold[hi], e_mm[hi], pch=17, col="#ff7f0e")
points(e_gold[hi], e_full[hi], pch=1, col="#d62728"); abline(0,1,lty=2,col="grey60")
legend("topleft",bty="n",cex=0.85,pch=c(19,17,1),col=c("#1f77b4","#ff7f0e","#d62728"),
       legend=c("RB-LOO","moment match","PSIS-LOO"))
dev.off(); cat("\nWrote figure: ereal_epilepsy_olre.png\n")

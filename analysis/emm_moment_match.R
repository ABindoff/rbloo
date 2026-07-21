# =====================================================================
# E-MM : moment-matching head-to-head (closes gate G2).
# Logistic random-intercept GLMM via brms (cmdstanr). On the high-k folds,
# compare vs reloo-exact gold:  PSIS-LOO  |  loo_moment_match  |  RB-LOO (quadrature).
# Pre-registered H5: RB >= MM accuracy at lower cost (RB 0 optimisations).
# =====================================================================
suppressMessages({library(brms); library(loo); library(posterior)})
options(mc.cores=4)

# ---- RB-LOO: marginalise b_j by 1-D quadrature given each base draw (alpha,sigma)
#      (identical estimator to gate_g2_glmm.R; the integrated-over-RE predictive) ----
rb_loo_glmm <- function(alpha, sig, y, g, J, gb=seq(-10,10,length.out=64)){
  S<-length(alpha); N<-length(y); L_rb<-matrix(0,S,N)
  for(i in 1:N){ j<-g[i]; oth<-setdiff(which(g==j),i)
    eta<-outer(alpha,gb,`+`); p<-plogis(eta)
    logp<- -0.5*outer(1/sig^2,gb^2) - outer(log(sig),rep(1,length(gb)))
    if(length(oth)) for(o in oth) logp<-logp + (if(y[o]==1) log(p) else log1p(-p))
    mx<-apply(logp,1,max); w<-exp(logp-mx); w<-w/rowSums(w)
    pred<-rowSums(w*(if(y[i]==1) p else 1-p)); L_rb[,i]<-log(pmax(pred,1e-300)) }
  lr<-suppressWarnings(loo(L_rb, r_eff=rep(1,N)))
  list(k=lr$diagnostics$pareto_k, elpd=lr$pointwise[,"elpd_loo"])
}
# structural leverage (GLM Fisher info) predictor, a-priori (C1)
struct_lev <- function(alpha, sig, bmat, y, g, J){
  pbar<-plogis(mean(alpha)+colMeans(bmat)[g]); info<-pbar*(1-pbar); su_h<-mean(sig)
  G_FFj<-sapply(1:J,function(j) sum(info[g==j])+1/su_h^2); info/G_FFj[g]
}
sp<-function(a,b) suppressWarnings(cor(a,b,method="spearman"))
auc<-function(s,l){ if(sum(l)==0||sum(l)==length(l)) return(NA); r<-rank(s)
  (sum(r[l==1])-sum(l==1)*(sum(l==1)+1)/2)/(sum(l==1)*sum(l==0)) }
rmse<-function(a,b) sqrt(mean((a-b)^2))

gen_data <- function(seed, J=60){
  set.seed(seed); nj<-sample(c(1,1,1,2,2,3),J,replace=TRUE); N<-sum(nj); g<-rep(1:J,nj)
  b<-rnorm(J,0,2.0); y<-rbinom(N,1,plogis(-0.2+b[g]))
  list(dat=data.frame(y=y,g=factor(g)), y=y, g=g, J=J, N=N, nsing=sum(nj==1))
}

extract_base <- function(fit, J){
  dm<-as_draws_df(fit)
  alpha<-as.numeric(dm[["b_Intercept"]]); sig<-as.numeric(dm[["sd_g__Intercept"]])
  bcols<-paste0("r_g[",1:J,",Intercept]"); bmat<-as.matrix(dm[,bcols]); colnames(bmat)<-NULL
  list(alpha=alpha, sig=sig, bmat=bmat)
}

REPS <- as.integer(Sys.getenv("EMM_REPS", "6"))
res_rows <- list(); fold_store <- list(); fit0 <- NULL

for(r in 1:REPS){
  d <- gen_data(seed=5+r, J=60)
  cat(sprintf("\n=== rep %d/%d : J=%d N=%d singletons=%d ===\n", r, REPS, d$J, d$N, d$nsing))

  if(is.null(fit0)){
    fit0 <- brm(y~1+(1|g), data=d$dat, family=bernoulli(),
                chains=4, iter=1000, refresh=0, seed=5+r,
                backend="rstan", save_pars=save_pars(all=TRUE),
                control=list(adapt_delta=0.95))
    fit <- fit0
  } else {
    fit <- update(fit0, newdata=d$dat, recompile=FALSE, refresh=0, seed=5+r)
  }

  ## PSIS-LOO
  lw <- suppressWarnings(loo(fit))
  k_full <- lw$diagnostics$pareto_k; e_full <- lw$pointwise[,"elpd_loo"]
  hi <- which(k_full>0.7); n_hi <- length(hi)
  cat(sprintf("  PSIS-LOO: #(k>0.7)=%d  max k=%.2f\n", n_hi, max(k_full)))
  if(n_hi==0){ cat("  no high-k folds; skipping rep for the head-to-head\n"); next }

  ## RB-LOO (0 optimisations)
  bs <- extract_base(fit, d$J)
  rb <- rb_loo_glmm(bs$alpha, bs$sig, d$y, d$g, d$J)
  h_struct <- struct_lev(bs$alpha, bs$sig, bs$bmat, d$y, d$g, d$J)

  ## Moment matching (SOTA cheaper-than-refit IS cure; k>0.7 obs get an optimisation)
  t_mm <- Sys.time()
  mm <- tryCatch(suppressWarnings(loo_moment_match(fit, loo=lw, k_threshold=0.7)),
                 error=function(e){ cat("  MM ERR:", conditionMessage(e), "\n"); NULL })
  mm_secs <- as.numeric(Sys.time()-t_mm, units="secs")
  if(!is.null(mm)){ e_mm<-mm$pointwise[,"elpd_loo"]; k_mm<-mm$diagnostics$pareto_k
    n_mm_opt <- n_hi   # MM performs a per-observation optimisation on each k>0.7 fold
    cat(sprintf("  MM: %d optimisations (%.0fs); #(k>0.7 after)=%d  max k=%.2f\n",
                n_mm_opt, mm_secs, sum(k_mm>0.7), max(k_mm)))
  } else { e_mm<-rep(NA,length(k_full)); k_mm<-rep(NA,length(k_full)); n_mm_opt<-NA }

  ## reloo gold : exact refit of every k>0.7 fold
  t_rl <- Sys.time()
  rl <- suppressWarnings(reloo(fit, loo=lw, k_threshold=0.7))
  rl_secs <- as.numeric(Sys.time()-t_rl, units="secs")
  e_gold <- rl$pointwise[,"elpd_loo"]
  cat(sprintf("  reloo: %d refits (%.0fs)\n", n_hi, rl_secs))

  ## per-fold accuracy on the high-k folds
  fold_store[[r]] <- data.frame(rep=r, obs=hi,
    e_gold=e_gold[hi], e_full=e_full[hi], e_mm=e_mm[hi], e_rb=rb$elpd[hi],
    k_full=k_full[hi], k_mm=k_mm[hi], k_rb=rb$k[hi])

  res_rows[[r]] <- data.frame(rep=r, N=d$N, n_hi=n_hi,
    C1_rho=sp(h_struct,k_full), C1_auc=auc(h_struct,as.numeric(k_full>0.7)),
    cure_rb=mean(rb$k[hi]<0.7), cure_mm=if(!is.null(mm)) mean(k_mm[hi]<0.7) else NA,
    rmse_psis=rmse(e_full[hi],e_gold[hi]),
    rmse_mm=if(!is.null(mm)) rmse(e_mm[hi],e_gold[hi]) else NA,
    rmse_rb=rmse(rb$elpd[hi],e_gold[hi]),
    opt_psis=0, opt_mm=n_mm_opt, opt_rb=0, refit_reloo=n_hi,
    mm_secs=mm_secs, reloo_secs=rl_secs)
  print(res_rows[[r]][,c("rep","n_hi","rmse_psis","rmse_mm","rmse_rb","opt_mm","opt_rb")])
}

RES <- do.call(rbind, res_rows); FOLD <- do.call(rbind, fold_store)
mc <- function(x){ x<-x[is.finite(x)]; c(mean(x), sd(x)/sqrt(length(x))) }

cat("\n================ E-MM : moment-matching head-to-head ================\n")
cat(sprintf("reps with high-k folds: %d ; total high-k folds: %d\n", nrow(RES), nrow(FOLD)))
cat(sprintf("C1  (a-priori) Spearman=%.2f+/-%.2f  AUC=%.2f+/-%.2f\n",
            mc(RES$C1_rho)[1],mc(RES$C1_rho)[2],mc(RES$C1_auc)[1],mc(RES$C1_auc)[2]))
cat(sprintf("cure fraction (k<0.7 after): RB=%.0f%%  MM=%.0f%%\n",
            100*mc(RES$cure_rb)[1], 100*mc(RES$cure_mm)[1]))
cat("\n--- accuracy vs reloo-exact gold on high-k folds (elpd RMSE, pooled) ---\n")
cat(sprintf("  PSIS-LOO : %.4f\n", rmse(FOLD$e_full, FOLD$e_gold)))
cat(sprintf("  MM       : %.4f\n", rmse(FOLD$e_mm[is.finite(FOLD$e_mm)], FOLD$e_gold[is.finite(FOLD$e_mm)])))
cat(sprintf("  RB-LOO   : %.4f\n", rmse(FOLD$e_rb, FOLD$e_gold)))
cat("\n--- per-rep RMSE (mean +/- MCSE) ---\n")
cat(sprintf("  PSIS=%.4f+/-%.4f  MM=%.4f+/-%.4f  RB=%.4f+/-%.4f\n",
            mc(RES$rmse_psis)[1],mc(RES$rmse_psis)[2],mc(RES$rmse_mm)[1],mc(RES$rmse_mm)[2],
            mc(RES$rmse_rb)[1],mc(RES$rmse_rb)[2]))
cat("\n--- COST (per rep) ---\n")
cat(sprintf("  RB optimisations=0  |  MM optimisations=%.0f (%.0fs)  |  reloo refits=%.0f (%.0fs)\n",
            mc(RES$opt_mm)[1], mc(RES$mm_secs)[1], mc(RES$refit_reloo)[1], mc(RES$reloo_secs)[1]))

saveRDS(list(RES=RES, FOLD=FOLD), "emm_moment_match.rds")

## ---- figure : the head-to-head ----
png("emm_moment_match.png", width=1600, height=540, res=135)
par(mfrow=c(1,3), mar=c(4.3,4.4,3.2,1))
yl<-range(c(FOLD$e_gold, FOLD$e_full, FOLD$e_rb, FOLD$e_mm), na.rm=TRUE)
plot(FOLD$e_gold, FOLD$e_rb, pch=19, col="#1f77b4", xlim=yl, ylim=yl, cex=0.8,
     xlab="reloo-exact elpd (gold)", ylab="LOO elpd", main="Accuracy on high-k folds")
points(FOLD$e_gold, FOLD$e_mm, pch=17, col="#ff7f0e", cex=0.8)
points(FOLD$e_gold, FOLD$e_full, pch=1, col="#d62728", cex=0.8)
abline(0,1,lty=2,col="grey55")
legend("topleft",bty="n",cex=0.85,pch=c(19,17,1),col=c("#1f77b4","#ff7f0e","#d62728"),
       legend=c("RB-LOO","moment match","PSIS-LOO"))
# RMSE bar
rm_ps<-rmse(FOLD$e_full,FOLD$e_gold); rm_mm<-rmse(FOLD$e_mm[is.finite(FOLD$e_mm)],FOLD$e_gold[is.finite(FOLD$e_mm)]); rm_rb<-rmse(FOLD$e_rb,FOLD$e_gold)
bp<-barplot(c(PSIS=rm_ps,MM=rm_mm,RB=rm_rb), col=c("#d62728","#ff7f0e","#1f77b4"),
        ylab="elpd RMSE vs exact refit", main="Accuracy (lower better)")
text(bp, c(rm_ps,rm_mm,rm_rb), sprintf("%.3f",c(rm_ps,rm_mm,rm_rb)), pos=3, xpd=NA, cex=0.9)
# cost
opt<-c(PSIS=0, MM=mean(RES$opt_mm,na.rm=TRUE), RB=0, reloo=mean(RES$refit_reloo))
bp2<-barplot(opt, col=c("#d62728","#ff7f0e","#1f77b4","grey40"),
        ylab="per-fold optimisations / refits", main="Cost (per replicate)")
text(bp2, opt, sprintf("%.0f",opt), pos=3, xpd=NA, cex=0.9)
dev.off(); cat("\nWrote figure: emm_moment_match.png\n")

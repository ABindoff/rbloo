# =====================================================================
# P2 + gate G2 : logistic random-intercept GLMM via rstanarm (precompiled).
# RB-LOO (analytic downdate + 1-D quadrature) vs PSIS-LOO vs reloo (exact gold).
# Tests C1 (leverage predicts k) and C2/H5 (RB cures + accuracy + cost) off Gaussian.
# =====================================================================
suppressMessages({library(rstanarm); library(loo)})
options(mc.cores=4); set.seed(5)

J<-60; nj<-sample(c(1,1,1,2,2,3),J,replace=TRUE); N<-sum(nj); g<-rep(1:J,nj)
b_true<-rnorm(J,0,2.0); y<-rbinom(N,1,plogis(-0.2+b_true[g]))
dat<-data.frame(y=y, g=factor(g))
cat(sprintf("GLMM: J=%d groups, N=%d obs, %d singletons\n", J,N,sum(nj==1)))

fit<-stan_glmer(y~1+(1|g), data=dat, family=binomial, chains=4, iter=1000,
                refresh=0, seed=5, adapt_delta=0.95)
dm<-as.matrix(fit); alpha<-dm[,"(Intercept)"]; sig<-sqrt(dm[,"Sigma[g:(Intercept),(Intercept)]"])
B<-dm[,paste0("b[(Intercept) g:",1:J,"]")]; S<-length(alpha)

## (B1) PSIS-LOO
lw<-loo(fit); k_full<-lw$diagnostics$pareto_k; e_full<-lw$pointwise[,"elpd_loo"]
cat(sprintf("PSIS-LOO: max k=%.2f  #(k>0.7)=%d / %d\n", max(k_full), sum(k_full>0.7), N))

## our RB-LOO: marginalise b_j by 1-D quadrature given each base draw (alpha,sigma)
gb<-seq(-10,10,length.out=64)
L_rb<-matrix(0,S,N)
for(i in 1:N){ j<-g[i]; oth<-setdiff(which(g==j),i)
  eta<-outer(alpha,gb,`+`); p<-plogis(eta)
  logp<- -0.5*outer(1/sig^2,gb^2) - outer(log(sig),rep(1,length(gb)))
  if(length(oth)) for(o in oth) logp<-logp + (if(y[o]==1) log(p) else log1p(-p))
  mx<-apply(logp,1,max); w<-exp(logp-mx); w<-w/rowSums(w)
  pred<-rowSums(w*(if(y[i]==1) p else 1-p)); L_rb[,i]<-log(pmax(pred,1e-300)) }
lr<-suppressWarnings(loo(L_rb,r_eff=rep(1,N))); k_rb<-lr$diagnostics$pareto_k; e_rb<-lr$pointwise[,"elpd_loo"]
cat(sprintf("RB-LOO:  max k=%.2f  #(k>0.7)=%d\n", max(k_rb), sum(k_rb>0.7)))

## (gold) reloo: exact refit of high-k folds via rstanarm k_threshold
t0<-Sys.time(); lg<-suppressWarnings(loo(fit, k_threshold=0.7)); reloo_secs<-as.numeric(Sys.time()-t0,units="secs")
e_gold<-lg$pointwise[,"elpd_loo"]; n_refit<-sum(k_full>0.7)

## C1: structural leverage (GLM Fisher info) predicts k_full
pbar<-plogis(mean(alpha)+colMeans(B)[g]); info<-pbar*(1-pbar); su_h<-mean(sig)
G_FFj<-sapply(1:J,function(j) sum(info[g==j])+1/su_h^2); h_struct<-info/G_FFj[g]
sp<-function(a,b) suppressWarnings(cor(a,b,method="spearman"))
auc<-function(s,l){ if(sum(l)==0||sum(l)==length(l)) return(NA); r<-rank(s); (sum(r[l==1])-sum(l==1)*(sum(l==1)+1)/2)/(sum(l==1)*sum(l==0)) }

hi<-which(k_full>0.7); rmse<-function(a,b) sqrt(mean((a-b)^2))
cat("\n================ GATE G2 (GLMM, rstanarm) ================\n")
cat(sprintf("C1  Spearman(structural leverage, k_full)=%.2f  AUC(k>0.7)=%.2f\n", sp(h_struct,k_full), auc(h_struct,as.numeric(k_full>0.7))))
cat(sprintf("C2  RB cures: PSIS #(k>0.7)=%d -> RB #(k>0.7)=%d   (mean k %.2f -> %.2f)\n", sum(k_full>0.7), sum(k_rb>0.7), mean(k_full), mean(k_rb)))
cat(sprintf("H5  accuracy vs reloo-gold on %d high-k folds (elpd RMSE):  PSIS-LOO=%.3f   RB-LOO=%.3f\n",
            length(hi), rmse(e_full[hi],e_gold[hi]), rmse(e_rb[hi],e_gold[hi])))
cat(sprintf("H5  COST:  reloo refits=%d models (%.0fs)  ;  RB-LOO refits=0 (quadrature over existing draws)\n", n_refit, reloo_secs))

saveRDS(list(k_full=k_full,k_rb=k_rb,e_full=e_full,e_rb=e_rb,e_gold=e_gold,h_struct=h_struct,hi=hi,y=y,g=g,
             alpha=alpha,sig=sig,B=B,dat=dat,fit_ok=TRUE), "gate_g2.rds")

png("gate_g2_glmm.png", width=1550, height=520, res=135)
par(mfrow=c(1,3), mar=c(4.2,4.3,3,1))
plot(h_struct,k_full,pch=19,col=adjustcolor("#333333",0.5),xlab="structural leverage (GLM Fisher, a-priori)",ylab="PSIS-LOO k-hat",
     main=sprintf("C1 on GLMM\nrho=%.2f AUC=%.2f",sp(h_struct,k_full),auc(h_struct,as.numeric(k_full>0.7))));abline(h=0.7,lty=3,col="red")
plot(k_full,k_rb,pch=19,col="#2ca02c",xlab="k-hat full (PSIS)",ylab="k-hat (RB-LOO)",main="C2: RB cures failures");abline(0,1,lty=2,col="grey60");abline(h=0.7,v=0.7,lty=3,col="grey70")
yl<-range(c(e_gold[hi],e_rb[hi],e_full[hi])); plot(e_gold[hi],e_rb[hi],pch=19,col="#1f77b4",xlim=yl,ylim=yl,xlab="reloo exact elpd (gold)",ylab="LOO elpd",main=sprintf("H5: accuracy (%d high-k folds)",length(hi)))
points(e_gold[hi],e_full[hi],pch=1,col="#d62728"); abline(0,1,lty=2,col="grey60"); legend("topleft",bty="n",cex=0.85,pch=c(19,1),col=c("#1f77b4","#d62728"),legend=c("RB-LOO","PSIS-LOO"))
dev.off(); cat("\nWrote figure: gate_g2_glmm.png\n")

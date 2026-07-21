# =====================================================================
# Replicate the logistic GLMM to put MCSE on C1/C2/H5. rstanarm (precompiled).
# (project: fibr_sbc)
# =====================================================================
suppressMessages({library(rstanarm); library(loo)})
options(mc.cores=4); set.seed(20)

rb_and_leverage <- function(fit, y, g, J){
  dm<-as.matrix(fit); alpha<-dm[,"(Intercept)"]; sig<-sqrt(dm[,"Sigma[g:(Intercept),(Intercept)]"])
  B<-dm[,paste0("b[(Intercept) g:",1:J,"]")]; S<-length(alpha); N<-length(y)
  lw<-suppressWarnings(loo(fit)); k_full<-lw$diagnostics$pareto_k; e_full<-lw$pointwise[,"elpd_loo"]
  gb<-seq(-10,10,length.out=64); L_rb<-matrix(0,S,N)
  for(i in 1:N){ j<-g[i]; oth<-setdiff(which(g==j),i)
    eta<-outer(alpha,gb,`+`); p<-plogis(eta)
    logp<- -0.5*outer(1/sig^2,gb^2)-outer(log(sig),rep(1,length(gb)))
    if(length(oth)) for(o in oth) logp<-logp+(if(y[o]==1) log(p) else log1p(-p))
    mx<-apply(logp,1,max); w<-exp(logp-mx); w<-w/rowSums(w)
    pred<-rowSums(w*(if(y[i]==1) p else 1-p)); L_rb[,i]<-log(pmax(pred,1e-300)) }
  lr<-suppressWarnings(loo(L_rb,r_eff=rep(1,N))); k_rb<-lr$diagnostics$pareto_k; e_rb<-lr$pointwise[,"elpd_loo"]
  pbar<-plogis(mean(alpha)+colMeans(B)[g]); info<-pbar*(1-pbar); su_h<-mean(sig)
  G_FFj<-sapply(1:J,function(j) sum(info[g==j])+1/su_h^2); h_struct<-info/G_FFj[g]
  list(h_struct=h_struct,k_full=k_full,k_rb=k_rb,e_full=e_full,e_rb=e_rb)
}
sp<-function(a,b) suppressWarnings(cor(a,b,method="spearman"))
auc<-function(s,l){ if(sum(l)==0||sum(l)==length(l)) return(NA); r<-rank(s); (sum(r[l==1])-sum(l==1)*(sum(l==1)+1)/2)/(sum(l==1)*sum(l==0)) }

REPS<-25; RELOO_REPS<-10
c1_rho<-c1_auc<-cured<-mkf<-mkr<-rep(NA,REPS)
poolfold<-list(); acc<-list()
for(r in 1:REPS){
  set.seed(100+r); J<-50; nj<-sample(c(1,1,1,2,2,3),J,replace=TRUE); N<-sum(nj); g<-rep(1:J,nj)
  b<-rnorm(J,0,2.0); y<-rbinom(N,1,plogis(-0.2+b[g])); dat<-data.frame(y=y,g=factor(g))
  fit<-suppressWarnings(stan_glmer(y~1+(1|g),data=dat,family=binomial,chains=4,iter=1000,refresh=0,seed=100+r,adapt_delta=0.95))
  o<-rb_and_leverage(fit,y,g,J)
  bad<-o$k_full>0.7
  c1_rho[r]<-sp(o$h_struct,o$k_full); c1_auc[r]<-auc(o$h_struct,as.numeric(bad))
  cured[r]<-if(any(bad)) mean(o$k_rb[bad]<0.7) else NA
  mkf[r]<-mean(o$k_full); mkr[r]<-mean(o$k_rb)
  poolfold[[r]]<-data.frame(h_struct=o$h_struct,k_full=o$k_full,k_rb=o$k_rb)
  if(r<=RELOO_REPS){ lw<-suppressWarnings(loo(fit)); lg<-suppressWarnings(loo(fit,k_threshold=0.7))
    eg<-lg$pointwise[,"elpd_loo"]; hi<-which(o$k_full>0.7)
    if(length(hi)) acc[[r]]<-data.frame(e_full=o$e_full[hi],e_rb=o$e_rb[hi],e_gold=eg[hi]) }
  cat(sprintf("  rep %2d/%d: k>0.7=%d  C1 rho=%.2f auc=%.2f  cured=%.0f%%\n",
              r,REPS,sum(bad),c1_rho[r],c1_auc[r],100*cured[r]))
}
PF<-do.call(rbind,poolfold); AC<-do.call(rbind,acc)
mc<-function(x) c(mean(x,na.rm=TRUE), sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))))
rmse<-function(a,b) sqrt(mean((a-b)^2))

cat("\n================ REPLICATED GLMM (mean +/- MCSE over replicates) ================\n")
r1<-mc(c1_rho); a1<-mc(c1_auc); cu<-mc(cured)
cat(sprintf("C1  Spearman(leverage,k_full) = %.3f +/- %.3f   |  AUC(k>0.7) = %.3f +/- %.3f\n", r1[1],r1[2],a1[1],a1[2]))
cat(sprintf("    pooled over %d folds: Spearman=%.3f  AUC=%.3f\n", nrow(PF), sp(PF$h_struct,PF$k_full), auc(PF$h_struct,as.numeric(PF$k_full>0.7))))
cat(sprintf("C2  cured fraction = %.1f%% +/- %.1f%%   (mean k_full %.2f -> k_rb %.2f)\n", 100*cu[1],100*cu[2], mean(mkf), mean(mkr)))
cat(sprintf("    pooled: %d folds with k_full>0.7 ; %d still >0.7 after RB\n", sum(PF$k_full>0.7), sum(PF$k_rb>0.7)))
cat(sprintf("H5  accuracy vs reloo-gold (%d reloo reps, %d high-k folds): RB RMSE=%.3f  PSIS RMSE=%.3f\n",
            RELOO_REPS, nrow(AC), rmse(AC$e_rb,AC$e_gold), rmse(AC$e_full,AC$e_gold)))
saveRDS(list(c1_rho=c1_rho,c1_auc=c1_auc,cured=cured,mkf=mkf,mkr=mkr,PF=PF,AC=AC),"replicate_glmm.rds")

png("replicate_glmm.png",width=1550,height=520,res=135)
par(mfrow=c(1,3),mar=c(4.2,4.3,3,1))
plot(PF$h_struct,PF$k_full,pch=19,col=adjustcolor("#333333",0.2),cex=0.5,xlab="structural leverage (a-priori)",ylab="PSIS-LOO k-hat",
     main=sprintf("C1 (pooled %d folds)\nAUC=%.2f",nrow(PF),auc(PF$h_struct,as.numeric(PF$k_full>0.7))));abline(h=0.7,lty=3,col="red")
plot(PF$k_full,PF$k_rb,pch=19,col=adjustcolor("#2ca02c",0.2),cex=0.5,xlab="k-hat full",ylab="k-hat RB",main="C2: RB cures (pooled)");abline(0,1,lty=2,col="grey60");abline(h=0.7,v=0.7,lty=3,col="grey70")
plot(AC$e_gold,AC$e_rb,pch=19,col="#1f77b4",cex=0.7,xlab="reloo exact elpd",ylab="LOO elpd",main=sprintf("H5 accuracy\nRB RMSE=%.3f PSIS=%.3f",rmse(AC$e_rb,AC$e_gold),rmse(AC$e_full,AC$e_gold)))
points(AC$e_gold,AC$e_full,pch=1,col="#d62728",cex=0.7);abline(0,1,lty=2,col="grey60");legend("topleft",bty="n",cex=0.85,pch=c(19,1),col=c("#1f77b4","#d62728"),legend=c("RB-LOO","PSIS-LOO"))
dev.off();cat("\nWrote figure: replicate_glmm.png\n")

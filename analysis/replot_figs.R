# Regenerate the trimmed figure set from saved .rds (no experiments re-run).
# Editorial pass: drop the repeated "RB cures" collapsed-cloud middle panels and
# the bar-charts-of-three-numbers; plain panel titles (no H-codes).
auc <- function(score,lab){ if(sum(lab)==0||sum(lab)==length(lab)) return(NA); r<-rank(score)
  (sum(r[lab==1])-sum(lab==1)*(sum(lab==1)+1)/2)/(sum(lab==1)*sum(lab==0)) }
rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))

## ---- Fig 1: Gaussian LMM (G1) -- two panels ----
g1 <- readRDS("rb_loo_gateG1.rds"); PF<-g1$PF; AC<-g1$AC
png("rb_loo_gateG1.png", width=1150, height=520, res=135)
par(mfrow=c(1,2), mar=c(4.3,4.4,2.6,1.4), cex.main=0.98)
a1 <- auc(PF$h_struct, as.numeric(PF$k_full>0.7))
plot(PF$h_struct, PF$k_full, pch=19, col=adjustcolor("#333333",0.25), cex=0.6,
     xlab="structural leverage (a priori)", ylab="PSIS-LOO k-hat",
     main=sprintf("Predicts k-hat (AUC %.2f)", a1)); abline(h=0.7,lty=3,col="#d62728")
plot(AC$e_ref, AC$e_rb, pch=19, col="#1f77b4", cex=0.7, xlab="exact refit elpd", ylab="LOO elpd",
     main=sprintf("vs refit (RMSE %.2f / %.2f)", rmse(AC$e_rb,AC$e_ref), rmse(AC$e_full,AC$e_ref)))
points(AC$e_ref, AC$e_full, pch=1, col="#d62728", cex=0.7); abline(0,1,lty=2,col="grey60")
legend("topleft",bty="n",cex=0.85,pch=c(19,1),col=c("#1f77b4","#d62728"),legend=c("RB-LOO","PSIS-LOO"))
dev.off(); cat("Fig G1 (2-panel) written\n")

## ---- Fig 3: E-MM -- single mechanism scatter (bars -> text) ----
mm <- readRDS("emm_moment_match.rds"); FO<-mm$FOLD
png("emm_moment_match.png", width=640, height=560, res=135)
par(mar=c(4.4,4.4,2.6,1))
yl<-range(c(FO$e_gold,FO$e_full,FO$e_rb,FO$e_mm),na.rm=TRUE)
plot(FO$e_gold, FO$e_rb, pch=19, col="#1f77b4", xlim=yl, ylim=yl, cex=0.8,
     xlab="exact refit elpd", ylab="LOO elpd", main="High-k folds vs exact refit")
points(FO$e_gold, FO$e_mm, pch=17, col="#ff7f0e", cex=0.8)
points(FO$e_gold, FO$e_full, pch=1, col="#d62728", cex=0.8); abline(0,1,lty=2,col="grey55")
legend("topleft",bty="n",cex=0.82,pch=c(19,17,1),col=c("#1f77b4","#ff7f0e","#d62728"),
       legend=c("RB-LOO","moment match","PSIS-LOO"))
dev.off(); cat("Fig E-MM (1-panel) written\n")

## ---- Fig 5: OLRE stressor -- two panels (drop cure cloud) ----
ol <- readRDS("ereal_epilepsy_olre.rds"); hi<-ol$hi
png("ereal_epilepsy_olre.png", width=1150, height=520, res=135)
par(mfrow=c(1,2), mar=c(4.3,4.4,2.6,1.4), cex.main=0.98)
a5 <- auc(ol$h_struct, as.numeric(ol$k_full>0.7))
plot(ol$h_struct, ol$k_full, pch=19, col=adjustcolor("#333333",0.5), cex=0.6,
     xlab="structural leverage (a priori)", ylab="PSIS-LOO k-hat",
     main=sprintf("Predicts k-hat (AUC %.2f)", a5)); abline(h=0.7,lty=3,col="#d62728")
yl<-range(c(ol$e_gold[hi],ol$e_rb[hi],ol$e_full[hi],ol$e_mm[hi]),na.rm=TRUE)
plot(ol$e_gold[hi], ol$e_rb[hi], pch=19, col="#1f77b4", xlim=yl, ylim=yl, cex=0.7,
     xlab="exact refit elpd", ylab="LOO elpd",
     main=sprintf("Accuracy on %d high-k folds", length(hi)))
points(ol$e_gold[hi], ol$e_mm[hi], pch=17, col="#ff7f0e", cex=0.7)
points(ol$e_gold[hi], ol$e_full[hi], pch=1, col="#d62728", cex=0.7); abline(0,1,lty=2,col="grey60")
legend("topleft",bty="n",cex=0.82,pch=c(19,17,1),col=c("#1f77b4","#ff7f0e","#d62728"),
       legend=c("RB-LOO","moment match","PSIS-LOO"))
dev.off(); cat("Fig OLRE (2-panel) written\n")

## ---- Fig 7: E-basestress -- two panels (drop the RMSE bar) ----
PA <- readRDS("e_basestress.rds")
flag <- PA$k_base>0.7
png("e_basestress.png", width=1150, height=520, res=135)
par(mfrow=c(1,2), mar=c(4.3,4.4,2.6,1.4), cex.main=0.98)
plot(PA$k_base, abs(PA$elpd_rb-PA$e_ref), pch=19, col=adjustcolor("#d62728",0.35), cex=0.6,
     xlab="RB-LOO base k-hat", ylab="| RB-LOO elpd  -  exact refit |",
     main="Error tracks base k-hat"); abline(v=0.7,lty=3,col="grey50")
cols<-ifelse(flag,"#d62728","#1f77b4")
plot(PA$e_ref, PA$elpd_rb, pch=19, col=adjustcolor(cols,0.5), cex=0.6,
     xlab="exact refit elpd", ylab="RB-LOO elpd", main="Matches unflagged, diverges flagged")
abline(0,1,lty=2,col="grey60")
legend("topleft",bty="n",cex=0.8,pch=19,col=c("#1f77b4","#d62728"),
       legend=c("base k-hat <= 0.7","base k-hat > 0.7 (refit)"))
dev.off(); cat("Fig basestress (2-panel) written\n")
cat("done.\n")

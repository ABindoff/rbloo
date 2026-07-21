# Final E-decision figure: (left) the model-selection verdict ladder;
# (right) the per-fold validation that RB = exact refit while PSIS drifts with k.
z <- readRDS("e_decision.rds"); D <- readRDS("adjudicate.rds")

png("e_decision.png", width=1550, height=560, res=135)
layout(matrix(c(1,2),1,2), widths=c(1.15,1))

## left: verdict ladder
par(mar=c(4.8,10,3.4,1.5))
labs <- c("RB-LOO\n(= exact, 0 refits)","reloo, k>0.7\n(standard remedy)","PSIS-LOO\n(naive default)")
# e_decision.rds fields: r_rb / r_relo / r_psis, each c(De=Delta, SE=SE)
D3  <- c(z$r_rb["De"],   z$r_relo["De"], z$r_psis["De"])
SE3 <- c(z$r_rb["SE"],   z$r_relo["SE"], z$r_psis["SE"])
cols<- c("#1f77b4","#7f7f7f","#d62728")
xlim<- c(-6, max(D3+2*SE3)+4)
plot(D3, 1:3, pch=19, col=cols, cex=1.7, yaxt="n", ylab="", xlim=xlim, ylim=c(0.5,3.5),
     xlab="elpd difference  (OLRE-Poisson  -  neg-binomial)",
     main="Model choice under three estimators")
segments(D3-2*SE3, 1:3, D3+2*SE3, 1:3, col=cols, lwd=3.5)
abline(v=0, lty=2, col="grey55")
axis(2, at=1:3, labels=labs, las=1, cex.axis=0.82)
text(D3, 1:3+0.24, sprintf("%+.1f +/- %.1f  (z=%.1f)", D3, 2*SE3, D3/SE3), cex=0.72, col=cols)
verdict <- c("indistinguishable","significant","decisive")
text(D3, 1:3-0.26, verdict, cex=0.72, font=3, col=cols)

## right: per-fold bias vs k -- RB unbiased, PSIS over-optimistic (grows with k)
par(mar=c(4.8,4.6,3.4,1.2))
yl <- range(c(D$err_psis, D$err_rb))
plot(D$k, D$err_psis, pch=1, col="#d62728", cex=0.9, ylim=yl, xlab="PSIS-LOO k-hat",
     ylab="LOO elpd  -  exact refit  (per fold)",
     main=sprintf("RB tracks exact refit  (RMSE %.3f vs %.3f)",
                  sqrt(mean(D$err_rb^2)), sqrt(mean(D$err_psis^2))))
points(D$k, D$err_rb, pch=19, col="#1f77b4", cex=0.75)
abline(h=0, lty=2, col="grey55"); abline(v=0.7, lty=3, col="grey60")
# trend of PSIS bias
lines(lowess(D$k, D$err_psis, f=0.8), col="#d62728", lwd=2)
legend("topleft", bty="n", cex=0.82, pch=c(19,1), col=c("#1f77b4","#d62728"),
       legend=c("RB-LOO - exact","PSIS-LOO - exact"))
text(0.72, yl[1]+0.05*diff(yl), "reloo(0.7)\nrefits only\nthis side ->", cex=0.6, col="grey45", adj=0)
dev.off(); cat("Wrote figure: e_decision.png\n")

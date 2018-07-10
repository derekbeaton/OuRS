rm(list=ls())
## some comparisons & tests.

library(tictoc)
library(ours)
library(GSVD)
library(rrcov)
library(robustbase)
library(cellWise)
library(golubEsets)


data("philips")
data("hbk")
data("Golub_Merge")

  ## my simple version
golub.data <- t(exprs(Golub_Merge))
golub.data.log10 <- log10(abs(golub.data)) * sign(golub.data)
golub.data.log10[is.nan(golub.data.log10)] <- 0

  ## from the second-generation p-value paper (see: https://github.com/LucyMcGowan/sgpvalue/blob/master/figures/GolubLeukemiaGeneData.R)
exprsDat <- exprs(Golub_Merge)
N.tmp <- nrow(exprsDat)
NormalizedGeneDat <-  apply(exprsDat, 2, function(z) qnorm((rank(z)-0.5)/N.tmp))
mostExtreme <- which(abs(NormalizedGeneDat) == max(abs(NormalizedGeneDat)))[1]
NormalizedGeneDat <- NormalizedGeneDat[-mostExtreme,]
NormalizedGeneDat_t <- t(NormalizedGeneDat)
N <- nrow(NormalizedGeneDat)
LeukDat <- cbind(NormalizedGeneDat_t,pData(Golub_Merge)[,c('ALL.AML','BM.PB','T.B.cell','Gender','PS','Source')])
hold <- LeukDat[with(LeukDat, order(LeukDat[7129])),]
leukdata <- hold[,-(7129:7134)]



print("START")


print("RRCOV MCD PHILIPS")
rrcov.mcd.philips_tic <- tic()
rrcov.mcd.philips <- CovMcd(philips)
rrcov.mcd.philips_toc <- toc()

print("RRCOV MCD HBK")
rrcov.mcd.hbk_tic <- tic()
rrcov.mcd.hbk <- CovMcd(hbk)
rrcov.mcd.hbk_toc <- toc()

print("OURS MCD PHILIPS")
ours.mcd.philips_tic <- tic()
ours.mcd.philips <- cont.mcd(philips)
ours.mcd.philips_toc <- toc()

print("OURS MCD HBK")
ours.mcd.hbk_tic <- tic()
ours.mcd.hbk <- cont.mcd(hbk)
ours.mcd.hbk_toc <- toc()


print("RRCOV HUBERT PHILIPS")
rrcov.hubert.philips_tic <- tic()
rrcov.hubert.philips <- PcaHubert(philips)
rrcov.hubert.philips_toc <- toc()

print("RRCOV HUBERT HBK")
rrcov.hubert.hbk_tic <- tic()
rrcov.hubert.hbk <- PcaHubert(hbk)
rrcov.hubert.hbk_toc <- toc()

print("RRCOV HUBERT LEUK")
rrcov.hubert.leuk_tic <- tic()
rrcov.hubert.leuk <- PcaHubert(leukdata)
rrcov.hubert.leuk_toc <- toc()

print("RRCOV HUBERT GOLUB")
rrcov.hubert.golub_tic <- tic()
rrcov.hubert.golub <- PcaHubert(golub.data.log10)
rrcov.hubert.golub_toc <- toc()

print("OURS SH PHILIPS")
ours.sh.philips_tic <- tic()
ours.sh.philips <- split.half.pca(philips)
ours.sh.philips_toc <- toc()

print("OURS SH HBK")
ours.sh.hbk_tic <- tic()
ours.sh.hbk <- split.half.pca(hbk)
ours.sh.hbk_toc <- toc()

print("OURS SH LEUK")
ours.sh.leuk_tic <- tic()
ours.sh.leuk <- split.half.pca(leukdata)
ours.sh.leuk_toc <- toc()

print("OURS SH GOLUB")
ours.sh.golub_tic <- tic()
ours.sh.golub <- split.half.pca(golub.data.log10)
ours.sh.golub_toc <- toc()

print("END")





## aggregate all the philips dists here

## there are 4 components
diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),mean))

source('tests/sh.dist.battery.R')

philips.pred.fi.mean <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,mean)
philips.pred.fi.median <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,median)
philips.pred.fi.max <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,max)
philips.pred.fi.iqr <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,IQR)
philips.pred.fi.sd <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,sd)
philips.pred.fi.diff <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})

philips.pred.u.mean <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,mean)
philips.pred.u.median <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,median)
philips.pred.u.max <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,max)
philips.pred.u.iqr <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,IQR)
philips.pred.u.sd <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,sd)
philips.pred.u.diff <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})

min.shr.det <- which(ours.sh.philips$sh.dets == min(ours.sh.philips$sh.dets), arr.ind = T)
if (min.shr.det[2] == 1) {
  shr.min.det.sample <- ours.sh.philips$sh1.orders[min.shr.det[1],
                                                           ]
} else {
  shr.min.det.sample <- ours.sh.philips$sh2.orders[min.shr.det[1],
                                                           ]
}
philips.shr.c.step <- cont.c.step(data = philips, obs.order = shr.min.det.sample, center = T, max.iters = 100)
philips.shr.c.step.pca <- tolerance.svd(scale(philips[philips.shr.c.step$obs.order, ], scale = F))
philips.shr.c.step.pred.fi <- scale(philips, center = colMeans(philips[philips.shr.c.step$obs.order, ]), scale = F) %*% philips.shr.c.step.pca$v
philips.shr.c.step.pred.u <- scale(philips, center = colMeans(philips[philips.shr.c.step$obs.order, ]), scale = F) %*% philips.shr.c.step.pca$v %*% diag(1/philips.shr.c.step.pca$d)
philips.shr.c.step.score <- rowSums(philips.shr.c.step.pred.fi^2)
philips.shr.c.step.mahal <- rowSums(philips.shr.c.step.pred.u^2)



sh.dist.res <- sh.dist.battery(philips,ours.sh.philips,center=T,scale=F)


philips.dists <- cbind(rrcov.mcd.philips@raw.mah,rrcov.mcd.philips@mah,ours.mcd.philips$dists$rob.md,ours.mcd.philips$dists$md,rrcov.hubert.philips@sd,rrcov.hubert.philips@od)
  colnames(philips.dists) <- c("rrcov.mcd.raw.mah","rrcov.mcd.mah","ours.mcd.rob.md","MD","rrcov.hubert.sd","rrcov.hubert.od")
philips.sh.dists <- cbind(philips.pred.fi.mean,philips.pred.fi.max,philips.pred.fi.iqr,philips.pred.fi.diff,philips.pred.fi.median,philips.pred.fi.sd,
                          philips.pred.u.mean,philips.pred.u.max,philips.pred.u.iqr,philips.pred.u.diff,philips.pred.u.median,philips.pred.u.sd,
                          philips.shr.c.step.score,philips.shr.c.step.mahal)

all.philips.dists <- cbind(philips.dists,philips.sh.dists)

#corrplot(cor(all.philips.dists),method="number")

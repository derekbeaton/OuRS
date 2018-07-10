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


print("START")


print("RRCOV MCD PHILIPS")
rrcov.mcd.philips_tic <- tic()
rrcov.mcd.philips <- CovMcd(philips)
rrcov.mcd.philips_toc <- toc()

print("OURS MCD PHILIPS")
ours.mcd.philips_tic <- tic()
ours.mcd.philips <- cont.mcd(philips)
ours.mcd.philips_toc <- toc()


print("RRCOV HUBERT PHILIPS")
rrcov.hubert.philips_tic <- tic()
rrcov.hubert.philips <- PcaHubert(philips)
rrcov.hubert.philips_toc <- toc()

print("OURS SH PHILIPS")
ours.sh.philips_tic <- tic()
ours.sh.philips <- split.half.pca(philips)
ours.sh.philips_toc <- toc()

print("END")





## aggregate all the philips dists here

## there are 4 components
diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),mean))

source('tests/sh.dist.battery.R')

# philips.pred.fi.mean <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,mean)
# philips.pred.fi.median <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,median)
# philips.pred.fi.max <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,max)
# philips.pred.fi.iqr <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,IQR)
# philips.pred.fi.sd <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,sd)
# philips.pred.fi.diff <- apply(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})
#
# philips.pred.u.mean <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,mean)
# philips.pred.u.median <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,median)
# philips.pred.u.max <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,max)
# philips.pred.u.iqr <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,IQR)
# philips.pred.u.sd <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,sd)
# philips.pred.u.diff <- apply(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})
#
min.shr.det <- which(ours.sh.philips$sh.dets == min(ours.sh.philips$sh.dets), arr.ind = T)
if (min.shr.det[2] == 1) {
  shr.min.det.sample <- ours.sh.philips$sh1.orders[min.shr.det[1],
                                                   ]
} else {
  shr.min.det.sample <- ours.sh.philips$sh2.orders[min.shr.det[1],
                                                   ]
}
# philips.shr.c.step <- cont.c.step(data = philips, obs.order = shr.min.det.sample, center = T, max.iters = 100)
# philips.shr.c.step.pca <- tolerance.svd(scale(philips[philips.shr.c.step$obs.order, ], scale = F))
# philips.shr.c.step.pred.fi <- scale(philips, center = colMeans(philips[philips.shr.c.step$obs.order, ]), scale = F) %*% philips.shr.c.step.pca$v
# philips.shr.c.step.pred.u <- scale(philips, center = colMeans(philips[philips.shr.c.step$obs.order, ]), scale = F) %*% philips.shr.c.step.pca$v %*% diag(1/philips.shr.c.step.pca$d)
# philips.shr.c.step.score <- rowSums(philips.shr.c.step.pred.fi^2)
# philips.shr.c.step.mahal <- rowSums(philips.shr.c.step.pred.u^2)
#


sh.dist.res <- sh.dist.battery(philips,ours.sh.philips,center=T,scale=F)
sh.dist.res_high.var <- sh.dist.battery(philips,ours.sh.philips,center=T,scale=F,comps=1:4)
sh.dist.res_low.var <- sh.dist.battery(philips,ours.sh.philips,center=T,scale=F,comps=5:9)


## OK so here, let's use the number of components and the sample for the min det to compute a new V matrix & center.

small.philips <- expo.scale(philips[shr.min.det.sample,],center=T,scale=F)
svd.res <- tolerance.svd(small.philips,nu=4,nv=4)
rob.center <- attributes(small.philips)$`scaled:center`

t.scores <- expo.scale(philips,center=rob.center,scale=F) %*% svd.res$v
test.this <- t.scores %*% t(svd.res$v)
test.od <- rowSums((expo.scale(philips,center=rob.center,scale=F) - test.this)^2)
test.sd <- rowSums(t.scores^2)
test.md <- rowSums(scale(t.scores)^2)
test.md.2 <- rowSums((t.scores %*% diag(1/svd.res$d[1:4]))^2)

test.dists <- cbind(test.od,test.sd,test.md,test.md.2)
  colnames(test.dists) <- c("t.OD","t.SD","t.MD","t.MD2")


philips.dists <- cbind(rrcov.mcd.philips@raw.mah,rrcov.mcd.philips@mah,ours.mcd.philips$dists$rob.md,ours.mcd.philips$dists$md,rrcov.hubert.philips@sd,rrcov.hubert.philips@od)
colnames(philips.dists) <- c("rrcov.mcd.raw.mah","rrcov.mcd.mah","ours.mcd.rob.md","MD","rrcov.hubert.sd","rrcov.hubert.od")
# philips.sh.dists <- cbind(philips.pred.fi.mean,philips.pred.fi.max,philips.pred.fi.iqr,philips.pred.fi.diff,philips.pred.fi.median,philips.pred.fi.sd,
#                            philips.pred.u.mean,philips.pred.u.max,philips.pred.u.iqr,philips.pred.u.diff,philips.pred.u.median,philips.pred.u.sd,
#                            philips.shr.c.step.score,philips.shr.c.step.mahal)

# all.philips.dists <- cbind(philips.dists,philips.sh.dists)

  ### I want some form of a flipped OD...
all.philips.dists <- cbind(philips.dists,
                           sh.dist.res$dists,
                           #sh.dist.res_high.var$lr.sh.dists,
                           #sh.dist.res_high.var$lr.dists,
                           test.dists
                           #sh.dist.res_low.var$lr.sh.dists,
                           #sh.dist.res_low.var$lr.dists
                           )

corrplot(cor(all.philips.dists),method="number")

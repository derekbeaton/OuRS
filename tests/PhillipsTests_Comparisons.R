rm(list=ls())
## some comparisons & tests.

library(tictoc)
library(ours)
library(GSVD)
library(rrcov)
library(robustbase)
library(cellWise)
library(golubEsets)
library(corrplot)

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
full.svd.res <- tolerance.svd(expo.scale(philips,center=T,scale=F))
full.fi <- full.svd.res$u %*% diag(full.svd.res$d)



small.philips <- expo.scale(philips[shr.min.det.sample,],center=T,scale=F)
rob.center <- attributes(small.philips)$`scaled:center`

svd.res <- tolerance.svd(small.philips,nu=4,nv=4)
svd.res2 <- tolerance.svd(small.philips)

t.scores <- expo.scale(philips,center=rob.center,scale=F) %*% svd.res$v
t.scores2 <- expo.scale(philips,center=rob.center,scale=F) %*% svd.res2$v

test.this <- t.scores %*% t(svd.res$v)
test.this2 <- t.scores2 %*% t(svd.res2$v)

test.od <- rowSums((expo.scale(philips,center=rob.center,scale=F) - test.this)^2)
test.od2 <- rowSums((expo.scale(philips,center=rob.center,scale=F) - test.this2)^2)

test.sd <- rowSums(t.scores^2)
test.md <- rowSums(scale(t.scores)^2)
test.md.2 <- rowSums((t.scores %*% diag(1/svd.res$d[1:4]))^2)

test.dists <- cbind(test.od,test.sd,test.md,test.md.2,test.od2)
  colnames(test.dists) <- c("t.OD","t.SD","t.MD","t.MD2","t.OD2")

# from https://github.com/cran/rospca/blob/master/R/Robpca.R
  # XRc <- X - matrix(rep(Xh.svd$center, times=n), nrow=n, byrow=TRUE)
  #
  # Xtilde <- XRc %*% Xh.svd$loadings[,1:k] %*% t(Xh.svd$loadings[,1:k])
  # Rdiff <- XRc - Xtilde
  # odh <- apply(Rdiff, 1, vecnorm)

## so robust centered data - robust centered based on components.
  DAT <- expo.scale(philips,scale=F)
  ROB.DAT <- expo.scale(philips,center=rob.center,scale=F)
  SUB.PROJ.ROB.DAT <- ROB.DAT %*% svd.res$v %*% t(svd.res$v)
  PROJ.ROB.DAT <- ROB.DAT %*% svd.res2$v %*% t(svd.res2$v)
  SUB.PROJ.DAT <- DAT %*% svd.res$v %*% t(svd.res$v)
  PROJ.DAT <- DAT %*% svd.res2$v %*% t(svd.res2$v)

  #apply(DAT-ROB.DAT,1,vecnorm) / rowSums((DAT-ROB.DAT)^2)

  # # apply(DAT - ROB.DAT,1,vecnorm) # dumb
  # od.1 <- apply(DAT - SUB.PROJ.ROB.DAT,1,vecnorm)
  # # apply(DAT - PROJ.ROB.DAT,1,vecnorm) # dumb
  # od.2 <- apply(DAT - SUB.PROJ.DAT,1,vecnorm)
  # # apply(DAT - PROJ.DAT,1,vecnorm) # dumb

  od.3 <- apply(ROB.DAT - SUB.PROJ.ROB.DAT,1,vecnorm)
  # apply(ROB.DAT - PROJ.ROB.DAT,1,vecnorm) # dumb
  od.4 <- apply(ROB.DAT - SUB.PROJ.DAT,1,vecnorm)
  # apply(ROB.DAT - PROJ.DAT,1,vecnorm) #dumb

  # od.5 <- apply(SUB.PROJ.ROB.DAT - PROJ.ROB.DAT,1,vecnorm)
  # # apply(SUB.PROJ.ROB.DAT - SUB.PROJ.DAT,1,vecnorm) # dumb
  # od.6 <- apply(SUB.PROJ.ROB.DAT - PROJ.DAT,1,vecnorm)
  #
  # od.7 <- apply(PROJ.ROB.DAT - SUB.PROJ.DAT,1,vecnorm)
  # # apply(PROJ.ROB.DAT - PROJ.DAT,1,vecnorm) # dumb
  #
  # od.8 <- apply(SUB.PROJ.DAT - PROJ.DAT,1,vecnorm)
  #
  # test.ods <- cbind(od.1,od.2,od.3,od.4,od.5,od.6,od.7,od.8)


philips.dists <- cbind(rrcov.mcd.philips@raw.mah,rrcov.mcd.philips@mah,ours.mcd.philips$dists$rob.md,ours.mcd.philips$dists$md,rrcov.hubert.philips@sd,rrcov.hubert.philips@od)
colnames(philips.dists) <- c("rrcov.mcd.raw.mah","rrcov.mcd.mah","ours.mcd.rob.md","MD","rrcov.hubert.sd","rrcov.hubert.od")
# philips.sh.dists <- cbind(philips.pred.fi.mean,philips.pred.fi.max,philips.pred.fi.iqr,philips.pred.fi.diff,philips.pred.fi.median,philips.pred.fi.sd,
#                            philips.pred.u.mean,philips.pred.u.max,philips.pred.u.iqr,philips.pred.u.diff,philips.pred.u.median,philips.pred.u.sd,
#                            philips.shr.c.step.score,philips.shr.c.step.mahal)

# all.philips.dists <- cbind(philips.dists,philips.sh.dists)

## two additional distance options:
  ## infinity MD & score
  infinity.sd <- rowSums(apply(ours.sh.philips$pred.fi.array^2,c(1,2),max))
  infinity.md <- rowSums(apply(ours.sh.philips$pred.u.array^2,c(1,2),max))

  ## partial OD
    ## scores from robust center & vectors - standard scores
      ## fill in zeros where necessary.
  u.scores_filled <- t.scores_filled <- matrix(0,nrow(full.fi),ncol(full.fi))
  t.scores_filled[1:nrow(t.scores),1:ncol(t.scores)] <- t.scores
  # t.scores_filled[1:nrow(t.scores),1:ncol(t.scores)] <- t.scores

  fi.strange.component.od <- rowSums((full.fi - t.scores_filled)^2)
  fi.partial.component.od <- rowSums((full.fi[,1:ncol(t.scores)] - t.scores)^2)
  fi.remaining.component.od <- rowSums(full.fi[,5:9]^2)

  test.more.dists <- cbind(infinity.md,infinity.sd,fi.strange.component.od,fi.partial.component.od,fi.remaining.component.od)

  ## I want some form of a flipped OD...
    ## maybe an OD ratio? or inverse OD diff (sum of squared diff)?

  ## also I need to be able to use the score and MDs from the two subsets of components...
  #test.dat1 <- DAT %*% svd.res$v %*% t(svd.res$v)
  test.dat1 <- DAT %*% full.svd.res$v[,1:4] %*% t(full.svd.res$v[,1:4])
  test.dat2 <- full.svd.res$u[,1:4] %*% diag(full.svd.res$d[1:4]) %*% t(full.svd.res$v[,1:4])


  last.od <- apply(DAT - test.dat2,1,vecnorm)

#test.res <- cont.mcd(test.dat2 + matrix(colMeans(philips),nrow(DAT),ncol(DAT),byrow=T),collinearity.stop = F)

all.philips.dists <- cbind(philips.dists,
                           sh.dist.res$dists,
                           #sh.dist.res_high.var$lr.sh.dists,
                           #sh.dist.res_high.var$lr.dists,
                           test.dists,
                           last.od
                           #rowSums(full.svd.res$u[,1:4]^2),rowSums(full.svd.res$u[,5:9]^2),rowSums(full.svd.res$u[,2:9]^2),rowSums(full.svd.res$u[,8:9]^2),
                           #sh.dist.res$dists[,c("pred.fi.median")] + sh.dist.res$dists[,c("pred.fi.iqr")],
                           #sh.dist.res$dists[,c("pred.u.median")] + sh.dist.res$dists[,c("pred.u.iqr")]
                           #test.res$dists$rob.md,test.res$dists$md,test.res$dists$rob.chid
                           #sh.dist.res_low.var$lr.sh.dists,
                           #sh.dist.res_low.var$lr.dists
                           )


corrplot(cor(all.philips.dists,method="spearman"),method="number")



### my main issue: I cannot depend on identifying a subsample that is good
  ## I need to avoid that at almost any cost...
    ## so maybe a OD with no sub sample but the reproducible components?


# epPCA(all.philips.dists[,c("pred.fi.median","pred.u.median")])
# epPCA(all.philips.dists[,c("pred.fi.median","pred.u.median","pred.fi.iqr","pred.u.iqr")])


test.pca <- epPCA(cbind(all.philips.dists[,c("pred.fi.median","pred.u.median")],last.od),graphs=F)


### OK I think the distances to work with are pred.fi.med or iqr, and pred.u.med or iqr , plus last.od
  ## additionally I think we can use the full set of distances from resampling to find a natural cutoff


score.distrs <- sqrt(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum))
u.distrs <- sqrt(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum))

# score.distrs > quantile(c(score.distrs),.95)
# u.distrs > quantile(c(u.distrs),.95)

bottom.score <- length(c(score.distrs)) * .025
top.score <- length(c(score.distrs)) * .975
top.score2 <- length(c(score.distrs)) * .95
up.low.95 <- sort(c(score.distrs))[round(c(bottom.score,top.score))]
up.95 <- sort(c(score.distrs))[round(c(top.score2))]

boxplot(t(score.distrs[order(apply(score.distrs,1,median)),]))
abline(h=up.95,col="red",lty=1)
abline(h=up.low.95[1],col="purple",lty=2)
abline(h=up.low.95[2],col="purple",lty=2)

score.cuts <- apply(score.distrs,1,function(x){sum(x < up.95) / length(x)})

bottom.u <- length(c(u.distrs)) * .025
top.u <- length(c(u.distrs)) * .975
top.u2 <- length(c(u.distrs)) * .95
u.up.low.95 <- sort(c(u.distrs))[round(c(bottom.u,top.u))]
u.up.95 <- sort(c(u.distrs))[round(c(top.u2))]

boxplot(t(u.distrs[order(apply(u.distrs,1,median)),]))
abline(h=u.up.95,col="red",lty=1)
abline(h=u.up.low.95[1],col="purple",lty=2)
abline(h=u.up.low.95[2],col="purple",lty=2)

u.cuts <- apply(u.distrs,1,function(x){sum(x < u.up.95) / length(x)})


  ## OK so there are two OD estimates
tolEllipsePlot(cbind(last.od,rowSums(tolerance.svd(S)$u^2)),classic=T)


  ## outliers can be based on predictions of U, FI, or how the two ODs come out.

# test.res <- apply(score.distrs,1,function(x){
#
#   x.sort <- sort(x)
#   bottom.x <- floor(length(x.sort)*.025)
#   top.x <- ceiling(length(x.sort)*.975)
#   sgpvalue::p_delta(x.sort[bottom.x],x.sort[top.x],0,up.95)
#
# })
#
# pca.test <- epPCA(cbind(apply(score.distrs,1,median),apply(u.distrs,1,median)),graphs=F)

# this.order <- order(pca.test$ExPosition.Data$fi[,1])

plot( cbind(apply(score.distrs,1,median),apply(u.distrs,1,median))[this.order,] )



  ## distribution-based outliers.
this.quant <- .95

score.outs <- rowSums(score.distrs > quantile(c(score.distrs), this.quant)) > (500 * this.quant)
m.outs <- rowSums(u.distrs > quantile(c(u.distrs), this.quant)) > (500 * this.quant)
od.outs <- last.od > quantile(last.od, this.quant)
fi.median.outs <- all.philips.dists[,c("pred.fi.median")] > quantile(all.philips.dists[,c("pred.fi.median")], this.quant)
u.median.outs <- all.philips.dists[,c("pred.u.median")] > quantile(all.philips.dists[,c("pred.u.median")], this.quant)


score.outs + m.outs + od.outs + fi.median.outs + u.median.outs

# (rowSums(score.distrs > quantile(c(score.distrs),.95)) > (500 * .95)) + (rowSums(u.distrs > quantile(c(u.distrs),.95)) > (500 * .95))
#
# which( ((rowSums(score.distrs > quantile(c(score.distrs),.95)) > (500 * .95)) + (rowSums(u.distrs > quantile(c(u.distrs),.95)) > (500 * .95))) > 0 )


#all.outs <- cbind( (rowSums(score.distrs > quantile(c(score.distrs),.95)) > (500 * .95)) , (rowSums(u.distrs > quantile(c(u.distrs),.95)) > (500 * .95)) , (last.od > quantile(last.od,.95)) )

#all.outs <- cbind(score.outs,m.outs,od.outs,fi.median.outs,u.median.outs)
all.outs <- cbind(score.outs,m.outs,od.outs)


intersect(which(!rrcov.hubert.philips@flag),which(rowSums(all.outs) > 0 ))
intersect(which(!rrcov.hubert.philips@flag),which(score.outs))
intersect(which(!rrcov.hubert.philips@flag),which(m.outs))
intersect(which(!rrcov.hubert.philips@flag),which(od.outs))
intersect(which(!rrcov.hubert.philips@flag),which(od.outs))
intersect(which(!rrcov.hubert.philips@flag),which(od.outs))



setdiff(which(!rrcov.hubert.philips@flag),which(rowSums(all.outs) > 0 ))


## OK I think I'm nearing a way to make decisions for cut-offs
  ## it should be kept simple, but I can point out that because of this framework, we have many other options

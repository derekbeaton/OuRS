
## Set up our pipeline to find outliers.
  ## first with the FI & U outliers.

  ## then as a post-hoc way, with robust subspace



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
library(ExPosition)

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

## ok so this is two different thresholds: one for where the line is, and one for the "probability of outlier-ness"
## both should be set to .95 by default, but the latter is not required to have a threshold.
## for visualization, the hard line can be drawn at the last point where there are non 0s
## then a second line at the cut point.
## this should be, effectively, a double boxplot with a shaded box.

dist.array.outliers <- function(dist.array, total.dist.cutoff = .95, outlier.cutoff = .95){

  dist.distrs <- sqrt(apply(dist.array^2,c(1,3),sum))
  upper.bound <- sort(c(dist.distrs))[round( length(c(dist.distrs)) * total.dist.cutoff )]

  outlier.scores <- apply(dist.distrs,1,function(x){sum(x >= upper.bound) / length(x)})
  outlier.threshold <- outlier.scores > outlier.cutoff

  return(list(dists=dist.distrs, cut.point = upper.bound, outlier.probabilities = outlier.scores, outliers = outlier.threshold))

}

sh.outliers <- function(sh.out, total.dist.cutoff = .95, outlier.cutoff = .95){

  score.outliers <- dist.array.outliers(sh.out$pred.fi.array)
  mahal.outliers <- dist.array.outliers(sh.out$pred.u.array)

  return(list(score.outliers=score.outliers,mahal.outliers=mahal.outliers))

}

reproducible.robust.low.rank.rebuild <- function(sh.out, corr.cutoff = NULL){

  # diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),mean))
  # diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),median))
  ## or an alternative would be to find out where the matrix smoothes out.

}

score.outlier.info <- dist.array.outliers(ours.sh.philips$pred.fi.array)
m.outlier.info <- dist.array.outliers(ours.sh.philips$pred.u.array)
sh.outlier.info <- sh.outliers(ours.sh.philips)


# median.dists <- cbind(apply(score.outlier.info$dists,1,median),apply(m.outlier.info$dists,1,median))
# max.dists <- cbind(apply(score.outlier.info$dists,1,max),apply(m.outlier.info$dists,1,max))
# min.dists <- cbind(apply(score.outlier.info$dists,1,min),apply(m.outlier.info$dists,1,min))

all.points <- cbind(c(sh.outlier.info$score.outliers$dists),c(sh.outlier.info$mahal.outliers$dists))
plot(0,0,type="n",xlim=c(0,max(all.points[,1])),ylim=c(0,max(all.points[,2])))
abline(v=sh.outlier.info$score.outliers$cut.point,col="olivedrab3")
abline(h=sh.outlier.info$mahal.outliers$cut.point,col="mediumorchid4")
for(i in 1:nrow(sh.outlier.info$score.outliers$dists)){
  this.ob <- cbind(sh.outlier.info$score.outliers$dists[i,],sh.outlier.info$mahal.outliers$dists[i,])

  if(sh.outlier.info$score.outliers$outliers[i]){
    if(sh.outlier.info$mahal.outliers$outliers[i]){
      this.col <- "firebrick3"
    }else{
      this.col <- "olivedrab3"
    }
  }else if(sh.outlier.info$mahal.outliers$outliers[i]){
    this.col <- "mediumorchid4"
  }else{
    this.col <- "grey80"
  }
  peeledHull(this.ob, col=this.col)
}

#
#
# obs.order <- order(tolerance.svd(expo.scale(median.dists))$u[,1])
#
# plot( median.dists[obs.order,] )
# abline(v = score.outlier.info$cut.point,col="red",lty=1,lwd=2)
# abline(h = m.outlier.info$cut.point,col="red",lty=1,lwd=2)

# ## first with the FI & U outliers.
#     ## THIS NEEDS A FUNCTION.
# score.distrs <- sqrt(apply(ours.sh.philips$pred.fi.array^2,c(1,3),sum))
# u.distrs <- sqrt(apply(ours.sh.philips$pred.u.array^2,c(1,3),sum))
#
# bottom.score <- length(c(score.distrs)) * .025
# top.score <- length(c(score.distrs)) * .975
# top.score2 <- length(c(score.distrs)) * .95
# up.low.95 <- sort(c(score.distrs))[round(c(bottom.score,top.score))]
# up.95 <- sort(c(score.distrs))[round(c(top.score2))]
#
# boxplot(t(score.distrs[order(apply(score.distrs,1,median)),]))
# abline(h=up.95,col="red",lty=1)
# abline(h=up.low.95[1],col="purple",lty=2)
# abline(h=up.low.95[2],col="purple",lty=2)
#
# score.cuts <- apply(score.distrs,1,function(x){sum(x >= up.95) / length(x)})
#
# bottom.u <- length(c(u.distrs)) * .025
# top.u <- length(c(u.distrs)) * .975
# top.u2 <- length(c(u.distrs)) * .95
# u.up.low.95 <- sort(c(u.distrs))[round(c(bottom.u,top.u))]
# u.up.95 <- sort(c(u.distrs))[round(c(top.u2))]
#
# boxplot(t(u.distrs[order(apply(u.distrs,1,median)),]))
# abline(h=u.up.95,col="red",lty=1)
# abline(h=u.up.low.95[1],col="purple",lty=2)
# abline(h=u.up.low.95[2],col="purple",lty=2)
#
# u.cuts <- apply(u.distrs,1,function(x){sum(x >= u.up.95) / length(x)})




## OK so there are two OD estimates
#diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),mean))
#diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),median))

mean.r2.mat <- apply(ours.sh.philips$loadings.cors^2,c(1,2),mean)
median.r2.mat <- apply(ours.sh.philips$loadings.cors^2,c(1,2),median)

block.r2 <- c()
for(i in 1:nrow(mean.r2.mat)){

  block.r2 <- c(block.r2,mean(c(mean.r2.mat[1:i,1:i])))

}


DAT <- expo.scale(philips,center=T,scale=F)
full.svd.res <- tolerance.svd(DAT)
low.rank.rebuild <- full.svd.res$u[,1:4] %*% diag(full.svd.res$d[1:4]) %*% t(full.svd.res$v[,1:4])

s.mat <- DAT - low.rank.rebuild
s.mat.svd <- tolerance.svd(s.mat)

s.mat.mds <- rowSums(s.mat.svd$u^2)
od <- apply(s.mat,1,vecnorm)


my.tol.plot <- function(dat,ellipse.alpha=.9,mcd.alpha=.9,graphs=F){

  pointsToEllipsoid <- function (X, Sigma, mu)
  {
    if (ncol(Sigma) != nrow(Sigma))
      stop("Sigma must be a square matrix")
    if (ncol(X) != ncol(Sigma))
      stop("number of columns in X must \n                                  be of same dimension as Sigma")
    if (length(mu) != ncol(Sigma))
      stop("length of mu must \n                                  be of same dimension as Sigma")
    eig <- eigen(Sigma)
    SigSqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    Z <- t(apply(X, 1, ellipsoidTransform, SigSqrt, mu))
    return(Z)
  }

  ## pTE private function. STOLEN FROM SIBER 2.1.0
  ellipsoidTransform <- function (x, SigSqrt, mu)
  {
    return(solve(SigSqrt, x - mu))
  }

  ## pTE private function. STOLEN FROM SIBER 2.1.0
  ellipseInOut <- function (Z, p = 0.95, r = NULL)
  {
    if (is.null(r)) {
      r <- stats::qchisq(p, df = ncol(Z))
    }
    inside <- rowSums(Z^2) < r
    return(inside)
  }

  addEllipse <- function (mu, sigma, m = NULL, n = 100, p.interval = NULL, ci.mean = FALSE,small.sample = FALSE, do.plot = TRUE, ...)
  {
    if (small.sample & is.null(m))
      message("A sample size number given by m is \n required when small.sample is TRUE")
    if (ci.mean & is.null(m))
      message("A sample size number given by m is \n  required when plotting confidence \n ellipses of the mean with ci.mean is TRUE")
    ifelse(ci.mean, c.scale <- m, c.scale <- 1)
    ifelse(small.sample, q <- (m - 1)/(m - 2), q <- 1)
    ifelse(is.null(p.interval), r <- 1, r <- sqrt(stats::qchisq(p.interval,
                                                                df = 2)))
    e = eigen(sigma/c.scale)
    SigSqrt = e$vectors %*% diag(sqrt(e$values * q)) %*% t(e$vectors)
    cc <- genCircle(n, r)
    back.trans <- function(x) {
      return(SigSqrt %*% x + mu)
    }
    ML.ellipse = t(apply(cc, 1, back.trans))
    if (grDevices::dev.cur() > 1 & do.plot) {
      graphics::lines(ML.ellipse, ...)
    }
    return(ML.ellipse)
  }

  genCircle <- function (n = 100, r)
  {
    theta = seq(0, 2 * pi, length = n)
    x = r * cos(theta)
    y = r * sin(theta)
    return(cbind(x, y))
  }

  #dat <- cbind(od,s.mat.mds)
  ## to get an estimate of just CD v MD robustness.
  mcd <- covMcd(dat,alpha = mcd.alpha)
  mcd.center <- mcd$center
  mcd.cov <- mcd$cov
  data.center <- colMeans(dat)
  data.cov <- cov(dat)

  z1 <- pointsToEllipsoid(dat,mcd.cov,mcd.center)
  z1.outs <- ellipseInOut(z1,p=ellipse.alpha)

  z2 <- pointsToEllipsoid(dat,data.cov,data.center)
  z2.outs <- ellipseInOut(z2,p=ellipse.alpha)


  if(graphs){
    x1 <- c(-max(dat[,1])*.05,max(dat[,1]))*1.1
    y1 <- c(-max(dat[,2])*.05,max(dat[,2]))*1.1


    plot(dat, xlim = x1, ylim = y1,pch=20,col="grey80", main=paste0("Outside of ellipse.alpha = ",ellipse.alpha), xlab="OD", ylab="Sub MD",cex=.5)
    tmp <- addEllipse(mcd.center,mcd.cov,p.interval = ellipse.alpha,col="red",lty=1)
    tmp2 <- addEllipse(data.center,data.cov,p.interval = ellipse.alpha,col="blue",lty=2)
    points(dat[!z1.outs,],bg="red",pch=21,cex=1)
    text(dat[!z1.outs,],labels=rownames(dat[!z1.outs,]),pos=1,col="red")
    points(dat[!z2.outs,],bg="blue",pch=21,cex=1)
    text(dat[!z2.outs,],labels=rownames(dat[!z2.outs,]),pos=1,col="blue")
    legend("bottomright",legend=c("Classic ellipse","Robust ellipse"), col=c("blue","red"), lty=c(2,1))
  }

  return(list(robust.outs=!z1.outs,classic.outs=!z2.outs))

}


#test.it <- cbind(ellipse.outliersb$robust.outs,ellipse.outliersb$classic.outs,sh.outlier.info$score.outliers$outliers,sh.outlier.info$mahal.outliers$outliers)+0

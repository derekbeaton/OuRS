
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
library(venneuler)
library(limma)



dist.array.outliers <- function(dist.array, total.dist.cutoff = .95, outlier.cutoff = .95){

  dist.array[is.na(dist.array)] <- 0

  dist.distrs <- sqrt(apply(dist.array^2,c(1,3),sum))
  upper.bound <- sort(c(dist.distrs))[round( length(c(dist.distrs)) * total.dist.cutoff )]

  outlier.scores <- apply(dist.distrs,1,function(x){sum(x >= upper.bound) / length(x)})
  outlier.threshold <- outlier.scores > outlier.cutoff

  return(list(dists=dist.distrs, cut.point = upper.bound, outlier.probabilities = outlier.scores, outliers = outlier.threshold))

}

sh.outliers <- function(sh.out, total.dist.cutoff = .95, outlier.cutoff = .95){

  score.outliers <- dist.array.outliers(sh.out$pred.fi.array, total.dist.cutoff = total.dist.cutoff, outlier.cutoff = outlier.cutoff)
  mahal.outliers <- dist.array.outliers(sh.out$pred.u.array, total.dist.cutoff = total.dist.cutoff, outlier.cutoff = outlier.cutoff)

  return(list(score.outliers=score.outliers,mahal.outliers=mahal.outliers))

}

reproducible.robust.low.rank.rebuild <- function(sh.out, corr.cutoff = NULL){

  # diag(apply(abs(ours.sh.leukdata$loadings.cors),c(1,2),mean))
  # diag(apply(abs(ours.sh.leukdata$loadings.cors),c(1,2),median))
  ## or an alternative would be to find out where the matrix smoothes out.

}

tol.ellipse <- function(dat,ellipse.alpha=.9,mcd.alpha=.9,xlab=colnames(dat)[1],ylab=colnames(dat)[2],graphs=F){

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

  genCircle <- function (n = 1000, r)
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

  rob.ellipse <- addEllipse(mcd.center,mcd.cov,p.interval = ellipse.alpha,col="blue",lty=2,do.plot = graphs)
  classic.ellipse <- addEllipse(data.center,data.cov,p.interval = ellipse.alpha,col="red",lty=2,do.plot = graphs)

  if(graphs){
    x1 <- c(-max(abs(dat[,1]))*.05,max(abs(dat[,1])))*1.1
    y1 <- c(-max(abs(dat[,2]))*.05,max(abs(dat[,2])))*1.1


    plot(dat, xlim = x1, ylim = y1,pch=20,col="grey80", main="", xlab=xlab, ylab=ylab,cex=.5)

    abline(v=max(rob.ellipse[,1]), h=max(rob.ellipse[,2]),lty=1,col="blue")
    points(dat[which(dat[,1] >= max(rob.ellipse[,1]) | dat[,2] >= max(rob.ellipse[,2])),],bg="blue",pch=21,cex=1)
    abline(v=max(classic.ellipse[,1]), h=max(classic.ellipse[,2]),lty=1,col="red")
    points(dat[which(dat[,1] >= max(classic.ellipse[,1]) | dat[,2] >= max(classic.ellipse[,2])),],bg="red",pch=21,cex=2)
    legend("bottomright",legend=c(paste0("Classic ellipse alpha = ", ellipse.alpha),paste0("Robust ellipse alpha = ", mcd.alpha)), col=c("red","blue"), lty=c(2,1))
  }

  return(
    list(
      x.robust.cutoff=max(rob.ellipse[,1]),
      x.classic.cutoff=max(classic.ellipse[,1]),
      y.robust.cutoff=max(rob.ellipse[,2]),
      y.classic.cutoff=max(classic.ellipse[,2]),

      x.robust.outliers = (dat[,1] >= max(rob.ellipse[,1])),
      x.classic.outliers= (dat[,1] >= max(classic.ellipse[,1])),
      y.robust.outliers = (dat[,2] >= max(rob.ellipse[,2])),
      y.classic.outliers= (dat[,2] >= max(classic.ellipse[,2]))
    )
  )

}


data("Golub_Merge")
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


print("RRCOV HUBERT leukdata")
rrcov.hubert.leukdata_tic <- tic()
rrcov.hubert.leukdata <- PcaHubert(leukdata)
rrcov.hubert.leukdata_toc <- toc()

print("OURS SH leukdata")
ours.sh.leukdata_tic <- tic()
ours.sh.leukdata <- split.half.pca(leukdata)
ours.sh.leukdata_toc <- toc()

print("END")

score.outlier.info <- dist.array.outliers(ours.sh.leukdata$pred.fi.array)
m.outlier.info <- dist.array.outliers(ours.sh.leukdata$pred.u.array)

score.outlier.info.sub <- dist.array.outliers(ours.sh.leukdata$pred.fi.array[,1:4,])
m.outlier.info.sub <- dist.array.outliers(ours.sh.leukdata$pred.u.array[,1:4,])

sh.outlier.info <- sh.outliers(ours.sh.leukdata)



## number of components requires inspection -- no way around it!
mean.r2.mat <- apply(ours.sh.leukdata$loadings.cors^2,c(1,2),mean)
median.r2.mat <- apply(ours.sh.leukdata$loadings.cors^2,c(1,2),median)
block.r2 <- c()
for(i in 1:nrow(mean.r2.mat)){

  block.r2 <- c(block.r2,mean(c(mean.r2.mat[1:i,1:i])))

}

#small.center <- colMeans(leukdata[ours.sh.leukdata$sh2.orders[96,],])
DAT <- expo.scale(leukdata,center=T,scale=F)
full.svd.res <- tolerance.svd(DAT)
low.rank.rebuild <- full.svd.res$u[,1:2] %*% diag(full.svd.res$d[1:2]) %*% t(full.svd.res$v[,1:2])

s.mat <- DAT - low.rank.rebuild
#s.mat.svd <- tolerance.svd(s.mat)

#s.mat.mds <- rowSums(s.mat.svd$u^2)
od <- apply(s.mat,1,vecnorm)





## distances to aggregate
  ## MD, robust MD, score d, orthogonal d; all from standard, MCD, ROBPCA
  ## median, IQR, 95% range for fi and u; od


  ## probably should sqrt these...
    ## all of these should actually come from sh.pca
      ## and I should throw in some component-wise ODs, as well.
      ## but for OD that can't work
      ## and I still think there needs to be this ellipse step.
my.dists <- cbind(
  apply(sqrt(score.outlier.info$dists),1,
        median),
  apply(sqrt(score.outlier.info$dists),1,
        IQR),
  apply(sqrt(score.outlier.info$dists),1,
        function(x){
          sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
        }),
  apply(sqrt(m.outlier.info$dists),1,
        median),
  apply(sqrt(m.outlier.info$dists),1,
        IQR),
  apply(sqrt(m.outlier.info$dists),1,
                       function(x){
                         sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
                       }),
  od
)

all.fin.dists <- cbind(rrcov.hubert.leukdata@od,rrcov.hubert.leukdata@sd,my.dists)
  colnames(all.fin.dists) <- c("ROBPCA OD","ROBPCA SD","SH SD median","SH SD IQR","SH SD 95%","SH MD median","SH MD IQR","SH MD 95%","SH OD")

  corrplot(cor(all.fin.dists),method="number")
  corrplot(cor(all.fin.dists,method = "spearman"),method="number")


  plot(rrcov.hubert.leukdata)
  plot(od, apply(m.outlier.info$dists,1,
                 function(x){
                   sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
                 }
  ))



md.interval.dists <- apply(m.outlier.info$dists,1,
                             function(x){
                               sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
                             }
        )

ellipse.data <- cbind(od,md.interval.dists)
  colnames(ellipse.data) <- c("od","md.intervals")



  te.res <- tol.ellipse(ellipse.data,graphs=T,ellipse.alpha = .5, mcd.alpha = .5)


    ## ok not bad... but also not great.
  all.three.method.outliers <- cbind(
    #(sqrt(rrcov.mcd.leukdata@raw.mah) >= mcd.cutoff)+0,
    (!rrcov.hubert.leukdata@flag)+0,
    (te.res$x.robust.outliers | te.res$y.robust.outliers)+0
  )

### I should be able to obtain the furthest point of the ellipse from 0... or just use quantiles?

  vennDiagram(vennCounts(all.three.method.outliers))


  #venn.diagram(list(ROBPCA=which((!rrcov.hubert.leukdata@flag)), SHPCA=which((te.res$x.robust.outliers | te.res$y.robust.outliers))),filename="test.tiff",fill=c("firebrick3","steelblue4"),alpha=c(.5,.5),cex=2,cat.fontface=4,lty=2,fontfamily=3)

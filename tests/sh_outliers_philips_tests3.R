
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
library(limma)

data("philips")


print("START")


plain.md <- mahalanobis(philips,colMeans(philips),cov(philips))

print("RRCOV MCD PHILIPS")
rrcov.mcd.philips_tic <- tic()
rrcov.mcd.philips <- CovMcd(philips)
rrcov.mcd.philips_toc <- toc()

print("RRCOV HUBERT PHILIPS")
rrcov.hubert.philips_tic <- tic()
rrcov.hubert.philips <- PcaHubert(philips)
rrcov.hubert.philips_toc <- toc()

print("OURS SH PHILIPS")
ours.sh.philips_tic <- tic()
ours.sh.philips <- split.half.pca(philips)
ours.sh.philips_toc <- toc()

print("OURS 75% PHILIPS")
ours.sh.philips75_tic <- tic()
ours.sh.philips75 <- two.fold.repeated.pca(philips,sh1.size = .75)
ours.sh.philips75_toc <- toc()


print("OURS 90% PHILIPS")
ours.sh.philips90_tic <- tic()
ours.sh.philips90 <- two.fold.repeated.pca(philips,sh1.size = .9)
ours.sh.philips90_toc <- toc()

print("END")


score.outlier.info_new <- make.distance.distributions.summaries(ours.sh.philips$pred.fi.array)
m.outlier.info_new <- make.distance.distributions.summaries(ours.sh.philips$pred.u.array)

score.outlier.info75_new <- make.distance.distributions.summaries(ours.sh.philips75$pred.fi.array)
m.outlier.info75_new <- make.distance.distributions.summaries(ours.sh.philips75$pred.u.array)

score.outlier.info90_new <- make.distance.distributions.summaries(ours.sh.philips90$pred.fi.array)
m.outlier.info90_new <- make.distance.distributions.summaries(ours.sh.philips90$pred.u.array)


## number of components requires inspection -- no way around it!
loadings.mean.r2.mat <- apply(ours.sh.philips$loadings.cors^2,c(1,2),mean)
loadings.median.r2.mat <- apply(ours.sh.philips$loadings.cors^2,c(1,2),median)

score.mean.r2.mat <- apply(ours.sh.philips$score.cors^2,c(1,2),mean)
score.median.r2.mat <- apply(ours.sh.philips$score.cors^2,c(1,2),median)

od_new <- low.rank.orthogonal.distances(philips,T,F,components=1:4)


my.dists2 <- cbind(
  score.outlier.info_new$median.dist,
  score.outlier.info_new$iqr.dist,
  score.outlier.info_new$percentile.dist,
  m.outlier.info_new$median.dist,
  m.outlier.info_new$iqr.dist,
  m.outlier.info_new$percentile.dist,
  od_new$od
)


score.outlier.scores <- sh.distribution.outliers(score.outlier.info_new$dists)
m.outlier.scores <- sh.distribution.outliers(m.outlier.info_new$dists)


all.fin.dists <- cbind(sqrt(plain.md),sqrt(rrcov.mcd.philips@raw.mah),sqrt(rrcov.mcd.philips@mah),rrcov.hubert.philips@od,rrcov.hubert.philips@sd,my.dists2)
  colnames(all.fin.dists) <- c("MD","MCD Robust MD","MCD Robust corrected MD","ROBPCA OD","ROBPCA SD","SH SD median","SH SD IQR","SH SD 95%","SH MD median","SH MD IQR","SH MD 95%","SH OD")

  corrplot(cor(all.fin.dists),method="number")
  corrplot(cor(all.fin.dists,method = "spearman"),method="number")


  plot(rrcov.mcd.philips)
  plot(rrcov.hubert.philips)


ellipse.data <- cbind(od_new$od,m.outlier.info_new$percentile.dist)
  colnames(ellipse.data) <- c("od","md.intervals")

### now also need a simple counting cutoff, like with the original dist outliers.
    ### just get X% of the distribution, and count how often each observation exists outside of that distribution.


  #### THIS IS THE WINNER FOR PRESENTATION.
      ### this actually does a fairly good job and exists somewhere in the middle.
      ### still need to note that we have other options.
  te.res <- tol.ellipse(ellipse.data,graphs=T)

  mcd.cutoff <- sqrt(qchisq(0.975, ncol(rrcov.mcd.philips@X)))
  all.outliers <- cbind(
    (sqrt(plain.md) >= mcd.cutoff)+0,
    (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)+0,
    (rrcov.hubert.philips@od >= rrcov.hubert.philips@cutoff.od)+0,
    (rrcov.hubert.philips@sd >= rrcov.hubert.philips@cutoff.sd)+0,
    (te.res$x.robust.outliers)+0,
    (te.res$y.robust.outliers)+0,
    (score.outlier.scores$outlier.scores > .95) + 0,
    (m.outlier.scores$outlier.scores > .95) + 0
  )
  crossprod(all.outliers)


  all.three.method.outliers <- cbind(
    (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)+0,
    (!rrcov.hubert.philips@flag)+0,
    (te.res$x.robust.outliers | te.res$y.robust.outliers)+0,
    (score.outlier.scores$outlier.scores > .95 | m.outlier.scores$outlier.scores > .95)+0
  )
  colnames(all.three.method.outliers) <- c("MCD outliers","ROBPCA outliers","SH PCA ellipse outliers","SH PCA distribution outliers")

  crossprod(all.three.method.outliers)
### I should be able to obtain the furthest point of the ellipse from 0... or just use quantiles?


# venn.diagram(list(which((sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)), which((!rrcov.hubert.philips@flag)), which((te.res$x.robust.outliers | te.res$y.robust.outliers))),filename = NULL)


  vennDiagram(vennCounts(all.three.method.outliers))


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
loadings75.median.r2.mat <- apply(ours.sh.philips75$loadings.cors^2,c(1,2),median)
loadings90.median.r2.mat <- apply(ours.sh.philips90$loadings.cors^2,c(1,2),median)

score.mean.r2.mat <- apply(ours.sh.philips$score.cors^2,c(1,2),mean)
score.median.r2.mat <- apply(ours.sh.philips$score.cors^2,c(1,2),median)
score.median75.r2.mat <- apply(ours.sh.philips75$score.cors^2,c(1,2),median)
score.median90.r2.mat <- apply(ours.sh.philips90$score.cors^2,c(1,2),median)

od_new <- low.rank.orthogonal.distances.test(philips,T,F,components=1:4, bootstrap.iters = 1000, alpha = .95, bootstrap.shortcut = F)


my.dists2 <- cbind(
  score.outlier.info_new$median.dist,
  score.outlier.info_new$iqr.dist,
  score.outlier.info_new$percentile.dist,
  m.outlier.info_new$median.dist,
  m.outlier.info_new$iqr.dist,
  m.outlier.info_new$percentile.dist,
  score.outlier.info75_new$median.dist,
  score.outlier.info75_new$iqr.dist,
  score.outlier.info75_new$percentile.dist,
  m.outlier.info75_new$median.dist,
  m.outlier.info75_new$iqr.dist,
  m.outlier.info75_new$percentile.dist,
  score.outlier.info90_new$median.dist,
  score.outlier.info90_new$iqr.dist,
  score.outlier.info90_new$percentile.dist,
  m.outlier.info90_new$median.dist,
  m.outlier.info90_new$iqr.dist,
  m.outlier.info90_new$percentile.dist,
  od_new$od
)


score.outlier.scores <- sh.distribution.outliers(score.outlier.info_new$dists)
m.outlier.scores <- sh.distribution.outliers(m.outlier.info_new$dists)

score.outlier.scores75 <- sh.distribution.outliers(score.outlier.info75_new$dists)
m.outlier.scores75 <- sh.distribution.outliers(m.outlier.info75_new$dists)

score.outlier.scores90 <- sh.distribution.outliers(score.outlier.info90_new$dists)
m.outlier.scores90 <- sh.distribution.outliers(m.outlier.info90_new$dists)


all.fin.dists <- cbind(sqrt(plain.md),sqrt(rrcov.mcd.philips@raw.mah),sqrt(rrcov.mcd.philips@mah),rrcov.hubert.philips@od,rrcov.hubert.philips@sd,my.dists2)

  colnames(all.fin.dists) <- c("MD","MCD Robust MD","MCD Robust corrected MD","ROBPCA OD","ROBPCA SD",paste0(rep(c("50","75","90"),each=6),c(" SH SD median","SH SD IQR","SH SD 95%","SH MD median","SH MD IQR","SH MD 95%")),"SH OD")

  corrplot(cor(all.fin.dists),method="number")
  corrplot(cor(all.fin.dists,method = "spearman"),method="number")

  # plot(rrcov.mcd.philips)
  # plot(rrcov.hubert.philips)

  mcd.cutoff <- sqrt(qchisq(0.975, ncol(rrcov.mcd.philips@X)))



  all.outliers <- cbind(
    (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)+0,
    (rrcov.hubert.philips@sd >= rrcov.hubert.philips@cutoff.sd)+0,
    (rrcov.hubert.philips@od >= rrcov.hubert.philips@cutoff.od)+0,
    (score.outlier.scores$outliers)+0,
    (m.outlier.scores$outliers)+0,
    (od_new$outliers)+0
  )
  colnames(all.three.method.outliers) <- c("MCD outliers","ROBPCA SD outliers","ROBPCA OD outliers","SH SD outliers","SH MD outliers","SH OD outliers")

  ### this is what we want/need. the distribution ones catch all the ROBPCA, substantial overlaps with MCD
      ### ok so both sets are important... we can catch all that exist in ROBPCA + a lot from MCD + others...
  all.three.method.outliers <- cbind(
    (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)+0,
    (!rrcov.hubert.philips@flag)+0,
    (score.outlier.scores$outliers  | m.outlier.scores$outliers | od_new$outliers)+0
  )
  colnames(all.three.method.outliers) <- c("MCD outliers","ROBPCA outliers","SH PCA distribution outliers")

  crossprod(all.three.method.outliers)
  vennDiagram(vennCounts(all.three.method.outliers))


#   ### ok so this is a good time to plot the MCD results and color them by (1) their outliers, (2) our outliers, and (3) overlap
#   pt.cols <- rep("grey80",nrow(philips))
#   pt.cols[(sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)] <- "olivedrab3"
#   pt.cols[(te.res$x.robust.outliers | te.res$y.robust.outliers | score.outlier.scores$outlier.scores > .95 | m.outlier.scores$outlier.scores > .95)] <- "mediumorchid4"
#   pt.cols[(te.res$x.robust.outliers | te.res$y.robust.outliers | score.outlier.scores$outlier.scores > .95 | m.outlier.scores$outlier.scores > .95) & (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)] <- "firebrick3"
#   plot(rrcov.mcd.philips,col=pt.cols,pch=20,labels="")
#     ## will need to re-do this plot.

  ## clearly the outlier cutoffs happen because of the 4
      ## I think I can simplify the ellipse approach but it is convenient and easy.

### I should be able to obtain the furthest point of the ellipse from 0... or just use quantiles?

  dist.mat <- cbind(od_new$od,score.outlier.info_new$median.dist,m.outlier.info_new$median.dist)
  colnames(dist.mat) <- c("SH OD","SH median SD","SH median MD")
  pca.res <- epPCA(dist.mat, DESIGN = (score.outlier.scores$outliers  | m.outlier.scores$outliers | od_new$outliers), make_design_nominal = T, graphs=F)
  epGraphs(pca.res,contributionPlots = F,correlationPlotter = F)


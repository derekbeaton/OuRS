
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


# dist.array.outliers <- function(dist.array, total.dist.cutoff = .95, outlier.cutoff = .95){
#
#   dist.distrs <- sqrt(apply(dist.array^2,c(1,3),sum))
#   upper.bound <- sort(c(dist.distrs))[round( length(c(dist.distrs)) * total.dist.cutoff )]
#
#   outlier.scores <- apply(dist.distrs,1,function(x){sum(x >= upper.bound) / length(x)})
#   outlier.threshold <- outlier.scores > outlier.cutoff
#
#   return(list(dists=dist.distrs, cut.point = upper.bound, outlier.probabilities = outlier.scores, outliers = outlier.threshold))
#
# }
#
# sh.outliers <- function(sh.out, total.dist.cutoff = .95, outlier.cutoff = .95){
#
#   score.outliers <- dist.array.outliers(sh.out$pred.fi.array, total.dist.cutoff = total.dist.cutoff, outlier.cutoff = outlier.cutoff)
#   mahal.outliers <- dist.array.outliers(sh.out$pred.u.array, total.dist.cutoff = total.dist.cutoff, outlier.cutoff = outlier.cutoff)
#
#   return(list(score.outliers=score.outliers,mahal.outliers=mahal.outliers))
#
# }
#
# reproducible.robust.low.rank.rebuild <- function(sh.out, corr.cutoff = NULL){
#
#   # diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),mean))
#   # diag(apply(abs(ours.sh.philips$loadings.cors),c(1,2),median))
#   ## or an alternative would be to find out where the matrix smoothes out.
#
# }




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

print("END")


score.outlier.info <- dist.array.outliers(ours.sh.philips$pred.fi.array)
m.outlier.info <- dist.array.outliers(ours.sh.philips$pred.u.array)


score.outlier.info_new <- make.distance.distributions.summaries(ours.sh.philips$pred.fi.array)
m.outlier.info_new <- make.distance.distributions.summaries(ours.sh.philips$pred.u.array)


## number of components requires inspection -- no way around it!
loadings.mean.r2.mat <- apply(ours.sh.philips$loadings.cors^2,c(1,2),mean)
loadings.median.r2.mat <- apply(ours.sh.philips$loadings.cors^2,c(1,2),median)

score.mean.r2.mat <- apply(ours.sh.philips$score.cors^2,c(1,2),mean)
score.median.r2.mat <- apply(ours.sh.philips$score.cors^2,c(1,2),median)



#small.center <- colMeans(philips[ours.sh.philips$sh2.orders[96,],])
DAT <- expo.scale(philips,center=T,scale=F)
full.svd.res <- tolerance.svd(DAT)
low.rank.rebuild <- full.svd.res$u[,1:4] %*% diag(full.svd.res$d[1:4]) %*% t(full.svd.res$v[,1:4])

s.mat <- DAT - low.rank.rebuild
s.mat.svd <- tolerance.svd(s.mat)

s.mat.mds <- rowSums(s.mat.svd$u^2)
od <- apply(s.mat,1,vecnorm)
my.od <- sqrt(rowSums(s.mat^2))

od_new <- low.rank.orthogonal.distances(philips,T,F,components=1:4)


## distances to aggregate
  ## MD, robust MD, score d, orthogonal d; all from standard, MCD, ROBPCA
  ## median, IQR, 95% range for fi and u; od


  ## probably should sqrt these...
my.dists <- cbind(
  apply(score.outlier.info$dists,1,
        median),
  apply(score.outlier.info$dists,1,
        IQR),
  apply(score.outlier.info$dists,1,
        function(x){
          sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
        }),
  apply(m.outlier.info$dists,1,
        median),
  apply(m.outlier.info$dists,1,
        IQR),
  apply(m.outlier.info$dists,1,
                       function(x){
                         sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
                       }),
  od
)

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
  # plot(od, apply(m.outlier.info$dists,1,
  #                function(x){
  #                  sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
  #                }
  # ))
  #
  # tol.ellipse.out <- tol.ellipse(cbind(od, apply(m.outlier.info$dists,1,
  #                      function(x){
  #                        sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
  #                      }
  # )),graphs = T,ellipse.alpha = .9, mcd.alpha = .9)
  #
  # tol.ellipse.out2 <- tol.ellipse(cbind(od, apply(score.outlier.info$dists,1,
  #                                            median
  # )),graphs = T)
  #
  # epPCA(cbind(od,apply(m.outlier.info$dists,1,
  #                      function(x){
  #                        sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
  #                      }
  # ),apply(score.outlier.info$dists,1,
  #         median
  # )))

#
#

# md.interval.dists <- apply(m.outlier.info$dists,1,
#                              function(x){
#                                sort(x)[ceiling(length(x)*.975)] - sort(x)[floor(length(x)*.025)]
#                              }
#         )

ellipse.data <- cbind(od_new$od,m.outlier.info_new$percentile.dist)
  colnames(ellipse.data) <- c("od","md.intervals")

### now also need a simple counting cutoff, like with the original dist outliers.
    ### just get X% of the distribution, and count how often each observation exists outside of that distribution.


  #### THIS IS THE WINNER FOR PRESENTATION.
      ### this actually does a fairly good job and exists somewhere in the middle.
      ### still need to note that we have other options.
  te.res <- tol.ellipse(ellipse.data,graphs=T,ellipse.alpha = .75, mcd.alpha = .75)

  mcd.cutoff <- sqrt(qchisq(0.975, ncol(rrcov.mcd.philips@X)))
  all.outliers <- cbind(
    (sqrt(plain.md) >= mcd.cutoff)+0,
    (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)+0,
    (rrcov.hubert.philips@od >= rrcov.hubert.philips@cutoff.od)+0,
    (rrcov.hubert.philips@sd >= rrcov.hubert.philips@cutoff.sd)+0,
    (te.res$x.robust.outliers)+0,
    (te.res$y.robust.outliers)+0
  )


  all.three.method.outliers <- cbind(
    (sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)+0,
    (!rrcov.hubert.philips@flag)+0,
    (te.res$x.robust.outliers | te.res$y.robust.outliers)+0
  )
  colnames(all.three.method.outliers) <- c("MCD outliers","ROBPCA outliers","SH PCA outliers")

  crossprod(all.three.method.outliers)
### I should be able to obtain the furthest point of the ellipse from 0... or just use quantiles?


# venn.diagram(list(which((sqrt(rrcov.mcd.philips@raw.mah) >= mcd.cutoff)), which((!rrcov.hubert.philips@flag)), which((te.res$x.robust.outliers | te.res$y.robust.outliers))),filename = NULL)

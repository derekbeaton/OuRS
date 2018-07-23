
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

score.outlier.info_new <- make.distance.distributions.summaries(ours.sh.leukdata$pred.fi.array)
m.outlier.info_new <- make.distance.distributions.summaries(ours.sh.leukdata$pred.u.array)


## number of components requires inspection -- no way around it!
loadings.mean.r2.mat <- apply(ours.sh.leukdata$loadings.cors^2,c(1,2),mean)
loadings.median.r2.mat <- apply(ours.sh.leukdata$loadings.cors^2,c(1,2),median)

score.mean.r2.mat <- apply(ours.sh.leukdata$score.cors^2,c(1,2),mean)
score.median.r2.mat <- apply(ours.sh.leukdata$score.cors^2,c(1,2),median)

#od_new <- low.rank.orthogonal.distances(leukdata,T,F,components=1:2)
od_new <- low.rank.orthogonal.distances.test(leukdata,T,F,components=1:2, bootstrap.iters = 1000, alpha = .95, bootstrap.shortcut = F)
#od_new_.75 <- low.rank.orthogonal.distances.test(leukdata,T,F,components=1:2, bootstrap.iters = 1000, alpha = .75, bootstrap.shortcut = F)


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

all.fin.dists <- cbind(rrcov.hubert.leukdata@od,rrcov.hubert.leukdata@sd,my.dists2)
  colnames(all.fin.dists) <- c("ROBPCA OD","ROBPCA SD","SH SD median","SH SD IQR","SH SD 95%","SH MD median","SH MD IQR","SH MD 95%","SH OD")

  corrplot(cor(all.fin.dists),method="number")
  corrplot(cor(all.fin.dists,method = "spearman"),method="number")

  plot(rrcov.hubert.leukdata)

  all.outliers <- cbind(
    (rrcov.hubert.leukdata@sd >= rrcov.hubert.leukdata@cutoff.sd)+0,
    (rrcov.hubert.leukdata@od >= rrcov.hubert.leukdata@cutoff.od)+0,
    (score.outlier.scores$outliers)+0,
    (m.outlier.scores$outliers)+0,
    (od_new$outliers)+0
  )
  colnames(all.outliers) <- c("ROBPCA SD outliers","ROBPCA OD outliers","SH SD outliers","SH MD outliers","SH OD outliers")

  ### this is what we want/need. the distribution ones catch all the ROBPCA, substantial overlaps with MCD
  ### ok so both sets are important... we can catch all that exist in ROBPCA + a lot from MCD + others...
  all.two.method.outliers <- cbind(
    (!rrcov.hubert.leukdata@flag)+0,
    (score.outlier.scores$outliers  | m.outlier.scores$outliers | od_new$outliers)+0
  )
  colnames(all.two.method.outliers) <- c("ROBPCA outliers","SH PCA distribution outliers")

  crossprod(all.two.method.outliers)
  vennDiagram(vennCounts(all.two.method.outliers))

    ## this isn't even really necessary...
  # dist.mat <- cbind(od_new$od,score.outlier.info_new$median.dist,m.outlier.info_new$median.dist)
  # colnames(dist.mat) <- c("SH OD","SH median SD","SH median MD")
  # pca.res <- epPCA(dist.mat, DESIGN = (score.outlier.scores$outliers  | m.outlier.scores$outliers | od_new$outliers), make_design_nominal = T, graphs=F)
  # epGraphs(pca.res,contributionPlots = F,correlationPlotter = F)

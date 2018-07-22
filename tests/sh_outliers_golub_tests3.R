
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

od_new <- low.rank.orthogonal.distances(leukdata,T,F,components=1:2)




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
  ellipse.data <- cbind(od_new$od,m.outlier.info_new$percentile.dist)
  colnames(ellipse.data) <- c("od","md.intervals")

  ### now also need a simple counting cutoff, like with the original dist outliers.
  ### just get X% of the distribution, and count how often each observation exists outside of that distribution.


  #### THIS IS THE WINNER FOR PRESENTATION.
  ### this actually does a fairly good job and exists somewhere in the middle.
  ### still need to note that we have other options.
  te.res <- tol.ellipse(ellipse.data,graphs=T)
  te.res_test <- tol.ellipse(sqrt(expo.scale(ellipse.data,scale="SS1")^2),graphs=T)


    ## ok not bad... but also not great.
  all.three.method.outliers <- cbind(
    #(sqrt(rrcov.mcd.leukdata@raw.mah) >= mcd.cutoff)+0,
    (!rrcov.hubert.leukdata@flag)+0,
    (te.res$x.robust.outliers | te.res$y.robust.outliers)+0,
    (score.outlier.scores$outlier.scores > .75 | m.outlier.scores$outlier.scores > .75)+0
  )
  colnames(all.three.method.outliers) <- c("ROBPCA outliers","SH PCA ellipse outliers","SH PCA distribution outliers")

### I should be able to obtain the furthest point of the ellipse from 0... or just use quantiles?

  vennDiagram(vennCounts(all.three.method.outliers))


  #venn.diagram(list(ROBPCA=which((!rrcov.hubert.leukdata@flag)), SHPCA=which((te.res$x.robust.outliers | te.res$y.robust.outliers))),filename="test.tiff",fill=c("firebrick3","steelblue4"),alpha=c(.5,.5),cex=2,cat.fontface=4,lty=2,fontfamily=3)


  all.two.method.outliers <- cbind(
    #(sqrt(rrcov.mcd.leukdata@raw.mah) >= mcd.cutoff)+0,
    (!rrcov.hubert.leukdata@flag)+0,
    (te.res$x.robust.outliers | te.res$y.robust.outliers | score.outlier.scores$outlier.scores > .75 | m.outlier.scores$outlier.scores > .75)+0
  )
  colnames(all.two.method.outliers) <- c("ROBPCA outliers","SH PCA  outliers")

  ### I should be able to obtain the furthest point of the ellipse from 0... or just use quantiles?

  vennDiagram(vennCounts(all.two.method.outliers))


  crossprod(cbind(te.res$x.robust.outliers+0, te.res$y.robust.outliers+0,(score.outlier.scores$outlier.scores > .75)+0 ,(m.outlier.scores$outlier.scores > .75)+0,(od_new$od > sort(sample(od_new$od,1000*length(od_new$od),replace=T))[ ceiling((1000*length(od_new$od))*.75) ])+0))






  low.rank.rebuild <- full.svd.res$u[,1:2] %*% diag(full.svd.res$d[1:2]) %*% t(full.svd.res$v[,1:2])

  DATA <- expo.scale(leukdata,center=T,scale=F)
  fixed.full.svd.res <- gsvd(DATA)

  all.ods <- matrix(0,nrow(leukdata),1000)
  for(i in 1:1000){
    boot.samp <- sample(nrow(leukdata),nrow(leukdata),replace=T)
    boot.data <- expo.scale(leukdata[boot.samp,],center=T,scale=F)

    low.rank.rebuild <- (boot.data %*% fixed.full.svd.res$v[,1:2]) %*% t(fixed.full.svd.res$v[,1:2])
    all.ods[,i] <- sqrt(rowSums( (boot.data - low.rank.rebuild) ^2))

  }



all.outs <-  cbind((od_new$od > sort(sample(od_new$od,1000*length(od_new$od),replace=T))[ ceiling((1000*length(od_new$od))*.75) ])+0,
        (od_new$od > sort(c(all.ods))[ceiling(length(c(all.ods)) * .75)])+0,
        te.res$x.robust.outliers+0,
        te.res$y.robust.outliers+0,
        #te.res$x.classic.outliers+0,
        #te.res$y.classic.outliers+0,
        (score.outlier.scores$outlier.scores > .75) + 0,
        (m.outlier.scores$outlier.scores > .75) + 0
        )



all.outs <-  cbind((!rrcov.hubert.leukdata@flag)+0,
                   ((od_new$od > sort(c(all.ods))[ceiling(length(c(all.ods)) * .75)]) |
                   (score.outlier.scores$outlier.scores > .75) + 0 |
                   (m.outlier.scores$outlier.scores > .75)) + 0
)

### ok so these are goig to be my final three...
  ## I need to add a bootstrap to low.rank.orthogonal, I suppose.


## so we use three criteria, inspired by ROBPCA & MCD
  ## but we can capture a proxy of MD.


### aggregate all three dists in a way into a PCA
  ## I guess OD as is, then maybe percentile
score.outlier.info_new2 <- make.distance.distributions.summaries(ours.sh.leukdata$pred.fi.array,lower.percentile = .125, upper.percentile = .875)
m.outlier.info_new2 <- make.distance.distributions.summaries(ours.sh.leukdata$pred.u.array,lower.percentile = .125, upper.percentile = .875)

  ## easy visual for this...
#epPCA(cbind(od_new$od,score.outlier.info_new2$median.dist,m.outlier.info_new2$median.dist), DESIGN = make.data.nominal(as.matrix(((od_new$od > sort(c(all.ods))[ceiling(length(c(all.ods)) * .75)]) |
                                                                                                                                            (score.outlier.scores$outlier.scores > .75) + 0 |
                                                                                                                                            (m.outlier.scores$outlier.scores > .75)))), make_design_nominal = F)

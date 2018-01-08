data(SNPS)

  ## will need these.
SNPs.nom <- make.data.nominal(SNPS)
ca.preproc.res <- ca.preproc(SNPs.nom)
profiles <- sweep(ca.preproc.res$Ox,1,ca.preproc.res$m,"/")

  ## do MCD
snps.res <- cat.mcd(SNPs.nom,make.data.disjunctive = F,num.subsets = 500,alpha = .75)

  ## do corrmax
snps.corrmax.res <- cat.corrmax(ca.preproc.res$weightedZx, loadings = snps.res$cov$loadings,singular.values = snps.res$cov$singular.values)

  ## do bootstrap
snps.boot.res <- cat.boot.sup.u(SNPs.nom,snps.res$cov$loadings,snps.res$cov$singular.values,iters = 100)


vec.boot.md <- c(snps.boot.res)

lower.cut <- sort(vec.boot.md)[(length(vec.boot.md) * .75)]
upper.cut <- sort(vec.boot.md)[(length(vec.boot.md) * .9)]

inliers <- which(snps.res$dists$rob.md <= lower.cut)
extreme.outliers <- which(snps.res$dists$rob.md >= upper.cut)
in.betweeniers <- which(snps.res$dists$rob.md < upper.cut & snps.res$dists$rob.md > lower.cut)


length(inliers) / nrow(SNPS)
length(extreme.outliers) / nrow(SNPS)
length(in.betweeniers) / nrow(SNPS)

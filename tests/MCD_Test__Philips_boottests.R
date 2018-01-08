library(cellWise)
data("philips")

ours.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)


## easy way out to identify problems:
## excluded from sample:
excluded <- setdiff(1:nrow(philips),ours.res$det.samps$samples[1,])



## OK so we do need a non-resampling cutoff... how can we do that? Just simple percentiles for now?

this.boot.samp <- 1:nrow(philips)
this.boot.samp <- c(1:676,676)

sup.scores <- cont.sup.fi.u(philips[this.boot.samp,],ours.res$cov$center,ours.res$cov$scale,loadings = ours.res$cov$loadings, singular.values = ours.res$cov$singular.values)

scale.dat <- expo.scale(philips[this.boot.samp,],T,F)
this.sup.u <- scale.dat %*% ours.res$cov$loadings %*% diag(1/ours.res$cov$singular.values)

  ## well, bootstrap doesn't really work here.
rowSums(this.sup.u^2) / ours.res$dists$rob.md[this.boot.samp]

test.set <- mvrnorm(nrow(philips), ours.res$cov$center, cov(philips[ours.res$det.samps$samples[1,],]))
#mvrnorm.mahal <- rowSums(tolerance.svd(test.set)$u^2)

mvnnorm.mahal.sup <- cont.sup.fi.u(test.set,T,F,ours.res$cov$loadings,ours.res$cov$singular.values)

#cbind(rowSums(this.sup.u^2),mvrnorm.mahal)


hist(rowSums(this.sup.u^2))
hist(mvrnorm.mahal)
hist(rowSums(mvnnorm.mahal.sup$sup.u^2))




### OK so I guess the bootstrap approach should not assume a robust center.

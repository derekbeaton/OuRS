data(SNPS)

snps.res <- cat.mcd(SNPS,make.data.disjunctive = T,num.subsets = 500,alpha = .5)
ca.preproc.res <- ca.preproc(make.data.nominal(SNPS))
profiles <- sweep(ca.preproc.res$Ox,1,ca.preproc.res$m,"/")

## OK so we do need a non-resampling cutoff... how can we do that? Just simple percentiles for now?

this.boot.samp <- 1:nrow(SNPS)
#this.boot.samp <- c(1:59,59)

#sup.scores <- cont.sup.fi.u(philips[this.boot.samp,],ours.res$cov$center,ours.res$cov$scale,loadings = ours.res$cov$loadings, singular.values = ours.res$cov$singular.values)
#sup.scores <- cat.sup.fi.u(profiles = profiles,row.weights = ca.preproc.res$m, col.weights = ca.preproc.res$w,loadings = ours.res$cov$loadings,singular.values = ours.res$cov$singular.values)

sup.res <- cat.sup.fi.u(profiles, ca.preproc.res$m, ca.preproc.res$w, snps.res$cov$loadings, snps.res$cov$singular.values)

scale.dat <- expo.scale(philips[this.boot.samp,],ours.res$cov$center,ours.res$cov$scale)
this.sup.u <- scale.dat %*% ours.res$cov$loadings %*% diag(1/ours.res$cov$singular.values)

## well, bootstrap doesn't really work here.
rowSums(this.sup.u^2) / ours.res$dists$rob.md[this.boot.samp]

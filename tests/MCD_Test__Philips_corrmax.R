library(cellWise)
data("philips")

ours.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)

corrmax.res <- cont.corrmax(target.data = philips,rob.center = ours.res$cov$center, rob.scale = ours.res$cov$scale,loadings = ours.res$cov$loadings, singular.values = ours.res$cov$singular.values)

heatmap(corrmax.res,Rowv=NA,Colv=NA)


## easy way out to identify problems:
  ## excluded from sample:
excluded <- setdiff(1:nrow(philips),ours.res$det.samps$samples[1,])


heatmap(corrmax.res[excluded,],Rowv=NA,Colv=NA)

  ## well that cutoff is quite clear.
plot(ours.res$dists$rob.md,ours.res$dists$u.od)


  ## OK so we do need a non-resampling cutoff... how can we do that? Just simple percentiles for now?

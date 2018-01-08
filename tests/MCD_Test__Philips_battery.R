library(cellWise)
data("philips")

  ## do MCD
philips.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)
      ### here we can use quantile or a diff cutoff.

  ## do corrmax -- use the robust center and scale
philips.corrmax.res <- cont.corrmax(target.data = philips,rob.center = philips.res$cov$center, rob.scale = philips.res$cov$scale,loadings = philips.res$cov$loadings, singular.values = philips.res$cov$singular.values)
      ### here we have percentages and can plot this as a heatmap.
        #### some sort of cutoffs should be used (e.g., 1/ncol(target.data) for variables; only outlier individuals)


  ## do bootstrap -- do NOT use the robust center and scale
philips.boot.res <- cont.boot.sup.u(philips,T,F,philips.res$cov$loadings,philips.res$cov$singular.values,iters = 100)
      ### here we have a boot distribution of possible distances, the matrix size is the same as the data but does not correspond to the rows
        #### this should be flattened/vectorized and used as larger distribution for a cutoff.



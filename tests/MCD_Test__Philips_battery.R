library(cellWise)
data("philips")

  ## do MCD
philips.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)

  ## do corrmax -- use the robust center and scale
philips.corrmax.res <- cont.corrmax(target.data = philips,rob.center = philips.res$cov$center, rob.scale = philips.res$cov$scale,loadings = philips.res$cov$loadings, singular.values = philips.res$cov$singular.values)

  ## do bootstrap -- do NOT use the robust center and scale
philips.boot.res <- cont.boot.sup.u(philips,T,F,philips.res$cov$loadings,philips.res$cov$singular.values,iters = 100)




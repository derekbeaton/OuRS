
library(robustbase)
#library(rrcov)
library(cellWise)
data("philips")

rb.res <- covMcd(philips,nsamp = 500,alpha = .75)
#rr.res <- CovMcd(philips,nsamp = 500)
ours.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)
  ## ours needs some speed up.

  ## we get the same sample.
#setdiff(rb.res$best,rr.res@best)
setdiff(rb.res$best,ours.res$best.sample)
setdiff(ours.res$best.sample,rb.res$best)
#setdiff(rr.res@best,ours.res$best.order)


  # what the hell is robustbase and rrcov doing to get Mahal?
#rb.res$raw.mah / rr.res@raw.mah
rb.res$raw.mah / ours.res$md
mahalanobis(philips,colMeans(philips),cov(philips)) / ours.res$md
mahalanobis(philips,colMeans(philips),cov(philips)) / rb.res$raw.mah


plot(rb.res$mah,rb.res$raw.mah) ## these are not "raw" Mahalanobis distances.
plot(rb.res$mah,ours.res$best.rob.md) ## we are not far off.
plot(rb.res$raw.mah,ours.res$best.rob.md) ## we are not far off.
plot(rb.res$raw.mah,ours.res$md) ## these are not "raw" Mahalanobis distances.

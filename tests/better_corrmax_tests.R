library(robustbase)
library(cellWise)
data("philips")


# better corrmax tests (the other script is a mess)

KS.corrmax <- function(target.data,rob.cov,rob.center){

  smpsd <- sqrt(diag(rob.cov))
  d <- diag(1/smpsd)      ## DB NOTE: a fraction can make a diag() panic.
  dsd <- d %*% rob.cov %*% d   ## DB Q: Is this converting cov to cor?

  trans1 <- sweep(target.data,2, rob.center)
  corrmax <- dsd %^% (-1/2) %*% d
  colnames(corrmax) <- rownames(corrmax) <- colnames(target.data)
  W <- t(apply(trans1,1,function(i){corrmax %*% as.matrix(i)}))
  #w2 <- t(contrib**2)

  return(sweep(W^2,1,rowSums(W^2),"/")*100)

}


rb.res <- covMcd(philips,nsamp = 500,alpha = .75)
ours.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)

setdiff(rb.res$best,ours.res$det.samps$samples[1,])
setdiff(ours.res$det.samps$samples[1,],rb.res$best)



RB_ks.corrmax.res <- KS.corrmax(target.data = philips,rob.cov = rb.res$cov,rob.center = rb.res$center)
RB_db.corrmax.res <- cont.corrmax(target.data = philips,rob.center = rb.res$center, rob.scale = F,loadings = eigen(rb.res$cov)$vectors, singular.values = sqrt(eigen(rb.res$cov)$values))

RB_db.corrmax.res / RB_ks.corrmax.res

OURS_ks.corrmax.res <- KS.corrmax(target.data = philips,rob.cov = ours.res$cov$loadings %*% diag(ours.res$cov$singular.values^2) %*% t(ours.res$cov$loadings),rob.center = ours.res$cov$center)
OURS_db.corrmax.res <- cont.corrmax(target.data = philips,rob.center = ours.res$cov$center, rob.scale = ours.res$cov$scale,loadings = ours.res$cov$loadings, singular.values = ours.res$cov$singular.values)

OURS_db.corrmax.res / OURS_ks.corrmax.res

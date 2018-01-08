library(cellWise)
data("philips")

ours.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)


## easy way out to identify problems:
## excluded from sample:
excluded <- setdiff(1:nrow(philips),ours.res$det.samps$samples[1,])



## OK so we do need a non-resampling cutoff... how can we do that? Just simple percentiles for now?

this.boot.samp <- 1:nrow(philips)
this.boot.samp <- c(1:676,676)
this.boot.samp <- sample(nrow(philips),replace=T)

sup.scores <- cont.sup.fi.u(philips[this.boot.samp,],ours.res$cov$center,ours.res$cov$scale,loadings = ours.res$cov$loadings, singular.values = ours.res$cov$singular.values)

  ##  guess the bootstrap part overall should be more about making decisions re: center/scale?
scale.dat <- expo.scale(philips[this.boot.samp,],T,F)
this.sup.u <- scale.dat %*% ours.res$cov$loadings %*% diag(1/ours.res$cov$singular.values)

  ## well, bootstrap doesn't really work here.
rowSums(this.sup.u^2) / ours.res$dists$rob.md[this.boot.samp]

#test.set <- mvrnorm(nrow(philips), ours.res$cov$center, cov(philips[ours.res$det.samps$samples[1,],]))
#mvrnorm.mahal <- rowSums(tolerance.svd(test.set)$u^2)

#mvnnorm.mahal.sup <- cont.sup.fi.u(test.set,T,F,ours.res$cov$loadings,ours.res$cov$singular.values)

#cbind(rowSums(this.sup.u^2),mvrnorm.mahal)


hist(rowSums(this.sup.u^2))
#hist(mvrnorm.mahal)
#hist(rowSums(mvnnorm.mahal.sup$sup.u^2))
hist(ours.res$dists$rob.md)




##real boot test
boot.maker <- function(target.data,center=T,scale=F,loadings,singular.values,iters=100){

  ### just go with the bootstrap approach but test what happens with very specific subsamples for the MD.
  boot.distrs <- matrix(NA,nrow(target.data),iters)
  for(i in 1:iters){

    boot.distrs[,i] <- rowSums( (expo.scale(target.data[sample(nrow(target.data),replace=T),],center=center,scale=scale) %*% loadings %*% diag(1/singular.values))^2 )

  }
  return(boot.distrs)
}




### OK so I guess the bootstrap approach should not assume a robust center.


#boot.res <- boot.maker(philips,T,F,ours.res$cov$loadings,ours.res$cov$singular.values,iters = 100)
boot.res <- cont.boot.sup.fi.u(philips,T,F,ours.res$cov$loadings,ours.res$cov$singular.values,iters = 100)

  ## compare these two.
hist(c(boot.res))
hist(ours.res$dists$rob.md)



#################
# BIG QUESTION
## What happens when we scale the SVs up by a constant?

## ANSWER:
  ## do not scale here. The scaling needs to be more of

### PERHAPS BIGGER ANSWER
  ## all it does, in our (OuRS; har har har) case is that it scales the MDs so no big deal
  ## however, the other scaling factors (in robustbase) that are more complex will have a different impact
#boot.res_constant <- boot.maker(philips,T,F,ours.res$cov$loadings,ours.res$cov$singular.values * (dim(ours.res$det.samps$samples)[2] / nrow(philips)),iters = 25)

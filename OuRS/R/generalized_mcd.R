
### the categorical version should explicitly handle the disjunctive transform
## same for ordinal and mixed data
## but all of these functions just pass the transformed data off to a generalized_mcd() which just takes data
  ## and it performs it as defined in the paper.


# ordinal_mcd <- function(){
#
# }


## data types will only allow for "cat" = categorical, "ord" for ordinal", "con" for continuous and "frq" for frequency
# mixed_data_mcd <- function(data, column.types=rep("cat",ncol(data)), alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){
# }


categorical_mcd <- function(data, make.data.disjunctive=F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  if(make.data.disjunctive){
    data <- make.data.nominal(data)
    make.data.disjunctive <- F
  }

  ### at this point it should call off to a generalized_mcd() which only takes data
    ## no passing of the disjunctive stuff.

  if(ncol(data) > (nrow(data)*.9)){
    stop("Column:Row ratio is too high. Please use another method (e.g., robPCA, rPCA, POET).")
  }
  if(alpha < .5){
    alpha <- .5
  }

  ## sample finder
  mcd.samples <- cat.mcd.find.sample(data, make.data.disjunctive=make.data.disjunctive, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent)
  ## only grab the top sample.
  best.sample <- mcd.samples$final.orders[1,]

  preproc.data <- ca.preproc(data) ## this could be more efficient...
  #profiles <- (diag(1/preproc.data$m) %*% preproc.data$Ox)
  profiles <- sweep(preproc.data$Ox,1,preproc.data$m,"/")

  ca.res <- ca(data)
  mahals <- rowSums(ca.res$u^2)
  chis <- rowSums(ca.res$fi^2)


  ## get robust mean & cov (loadings)
  rob.sample <- preproc.data$weightedZx[best.sample,]
  robust.tsvd.res <- tolerance.svd(rob.sample,tol=tol)

  ## call to function that computes robust mahal
  robust.dists <- cat.sup.fi.u(profiles, preproc.data$m, preproc.data$w, robust.tsvd.res$v, robust.tsvd.res$d)
  robust.mahals <- rowSums(robust.dists$sup.u^2)
  robust.chis <- rowSums(robust.dists$sup.fi^2)

  res <- list(
    cov = list(loadings = robust.tsvd.res$v,
               singular.values = robust.tsvd.res$d
    ),
    dists = list(rob.md = robust.mahals,
                 rob.chid = robust.chis,
                 md = mahals,
                 chid = chis),
    det.samps = list(dets = mcd.samples$final.dets,
                     samples = mcd.samples$final.orders)
  )


  class(res) <- c("list", "catMCD")
  return(res)

}


### techinically, the generalized MCD is CA-based, so it has to assume we have data that are CA-friendly
generalized_mcd <- function(data, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

    ## it's basically all of the above but without the transforms.


}

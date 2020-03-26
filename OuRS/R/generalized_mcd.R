
### the categorical version should explicitly handle the disjunctive transform
## same for ordinal and mixed data
## but all of these functions just pass the transformed data off to a generalized_mcd() which just takes data
  ## and it performs it as defined in the paper.


ordinal_mcd <- function(data, make_data_doubled=T, mins = NULL, maxs = NULL, impute_NA_to_mean = F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){


  if(make_data_doubled){
    data <- thermometer_coding(data, mins = mins, maxs = maxs, impute_NA_to_mean = impute_NA_to_mean)
  }

  generalized_mcd(data, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent, tol=tol)


}


## data types will only allow for "cat" = categorical, "ord" for ordinal", "con" for continuous and "frq" for frequency
# mixed_data_mcd <- function(data, column.types=rep("cat",ncol(data)), alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){
# }


categorical_mcd <- function(data, make_data_disjunctive=T, impute_NA_to_mean = F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  if(make_data_disjunctive){
    data <- disjunctive_coding(data, impute_NA_to_mean = impute_NA_to_mean)
  }

  generalized_mcd(data, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent, tol=tol)

}



### techinically, the generalized MCD is CA-based, so it has to assume we have data that are CA-friendly
generalized_mcd <- function(data, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

    ## it's basically all of the above but without the transforms.

  ### at this point it should call off to a generalized_mcd() which only takes data
  ## no passing of the disjunctive stuff.

  if(ncol(data) > (nrow(data)*.9)){
    stop("generalized_mcd: the column-to-row ratio is too high so 'mcd' cannot be performed")
  }
  if(alpha < .5){
    alpha <- .5
  }

  ## sample finder
  mcd.samples <- generalized_mcd_find_sample(data, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent)
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
  robust.dists <- generalized_scores_dists(profiles, preproc.data$m, preproc.data$w, robust.tsvd.res$v, robust.tsvd.res$d)
    ## change this
  # robust.mahals <- rowSums(robust.dists$sup.u^2)
  # robust.chis <- rowSums(robust.dists$sup.fi^2)

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
#### so this has the same requirement as generalized_mcd and works under the assumption that the data have been transformed appropriately

# formerly cat.mcd.find.sample
generalized_mcd_find_sample <- function(data, alpha=.75,num.subsets=500,max.total.iters=num.subsets*20,top.sets.percent=.05,tol=.Machine$double.eps){


  ## I need to assess how necesarry these various (possibly redundant) pieces are
  ## just compute these directly from data
  preproc.data <- ca.preproc(data) ## this could be more efficient...
  #profiles <- (diag(1/preproc.data$m) %*% preproc.data$Ox)
  profiles <- sweep(preproc.data$Ox,1,preproc.data$m,"/")

  if(alpha<.5){
    h.size <- floor((nrow(data)+1)/2)
  }
  h.size <- h.alpha.n(alpha,nrow(data),ncol(data))
  max.det.iters <- round(max.total.iters / num.subsets)

  dets <- vector("numeric", num.subsets)
  orders <- matrix(NA,num.subsets,h.size)
  for(i in 1:num.subsets){

    findInit <- T
    init.size <- min(dim(data))+1

    while( findInit ){

      init.samp <- sort(sample(nrow(data),init.size))

      init.svd <- tolerance.svd(preproc.data$weightedZx[init.samp,],tol=tol)
      init.mds <- round(rowSums(init.svd$u^2),digits=8)	## do I need to round? ### I should probably use a tol parameter here...

      if(length(unique(init.mds)) < 2){
        init.size <- init.size + 1
      }else{
        sup.scores <- generalized_scores_dists(profiles, preproc.data$m, preproc.data$w, init.svd$v, init.svd$d)
        ## change this
        # mahals <- rowSums(sup.scores$sup.u^2)
        samp.config <- sort(order(mahals)[1:h.size])
        findInit <- F
      }
    }

    min.info <- generalized_c_step(profiles, preproc.data$weightedZx, preproc.data$m, preproc.data$w, samp.config, max.det.iters)
    dets[i] <- min.info$min.det
    orders[i,] <- min.info$obs.order
  }

  perc.cut <- round(num.subsets * top.sets.percent)	## they take "top 10" -- I will allow our search to be broader
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique(unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),])

  final.dets <- vector("numeric", nrow(final.configs))
  final.orders <- matrix(NA,nrow(final.configs),h.size)
  for( i in 1:nrow(final.configs)){

    ## RETURN TO THIS. I NEED INF to work.
    min.info <- generalized_c_step(profiles, preproc.data$weightedZx, preproc.data$m, preproc.data$w, final.configs[i,], 1000)		## set to Inf so that this converges on its own; need to make this settable & have a real max embedded in c.step
    final.dets[i] <- min.info$min.det
    final.orders[i,] <- min.info$obs.order
  }

  best.order <- order(final.dets)
  final.dets <- final.dets[best.order]
  final.orders <- final.orders[best.order,]

  return( list(final.dets= final.dets, final.orders= final.orders) )

}

## same as the above two: this needs to be the generalized version, it doesn't care what data they were
  ## but also this should handle profiles internally
# formerly cat.sup.fi.u
generalized_scores_dists <- function(profiles, row.weights, col.weights, loadings, singular.values){

  ## NOT EFFICIENT. MAKE MORE EFFICIENT
  ## also this can handle the profiles here...

  sup.fi <- profiles %*% sweep(loadings,1,sqrt(col.weights)/col.weights,"*")
  sup.u <- sweep(sweep(sup.fi,2,singular.values,"/"),1,sqrt(row.weights),"*")

  mahals <- rowSums(sup.u^2)
  chis <- rowSums(sup.fi^2)

  return( list(mahals = mahals, chis = chis, sup.u=sup.u, sup.fi=sup.fi) )

}

### same same same as the above points: this is the generalized version

# formerly cat.c.step
generalized_c_step <- function(profiles, weighted.deviations, row.weights, col.weights, obs.order, max.iters=25, tol=sqrt(.Machine$double.eps)){

  old.det <- Inf
  old.center <- NaN
  old.v <- matrix(NaN,0,0)
  new.order <- old.order <- obs.order

  for(i in 1:max.iters){
    sub.data <- weighted.deviations[new.order,]
    new.center <- colMeans(sub.data)
    svd.res <- tryCatch( {tolerance.svd(sub.data)}, error=function(x) 'FAIL') ## this should probably be explained a bit!

    if(length(svd.res)==3){
      new.det <- geometric_mean(svd.res$d^2)


      if( (new.det <= old.det) & (!isTRUE(all.equal(new.det,0,tolerance=tol))) ){
        if( center.sigma_checker(old.center, new.center, old.v, svd.res$v,tol=tol) & isTRUE(all.equal(new.det, old.det, tolerance= tol)) ){
          return( list(obs.order = old.order, min.det = old.det) )
        }
        else{

          old.det <- new.det
          old.center <- new.center
          old.v <- svd.res$v
          old.order <- new.order

          sup.scores <- generalized_scores_dists(profiles, row.weights, col.weights, svd.res$v, svd.res$d)
          ## change this
          # mahals <- rowSums(sup.scores$sup.u^2)
          new.order <- sort(order(mahals)[1:length(obs.order)]	)

          if( isTRUE(all.equal(sort(new.order),sort(old.order))) ){
            return( list(obs.order = old.order, min.det = old.det) )
          }
        }
      }else{
        return( list(obs.order = old.order, min.det = old.det) )
      }
    }else{
      warning('SVD error caught and skipped. Prior results being returned.')
      return( list(obs.order = old.order, min.det = old.det) )
    }
  }
  ## this is basically the max out which implies the most recent is the best.
  return( list(obs.order = new.order, min.det= new.det) )
}



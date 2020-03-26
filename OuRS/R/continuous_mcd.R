
##
#'
#' @export
#'

continuous_mcd <- function(data, center=T, scale=F, collinearity.stop=T, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  ## data size check
  if(ncol(data) > (nrow(data)*.9)){
  	stop("continuous_mcd: the column-to-row ratio is too high so 'mcd' cannot be performed")
  }

  ## collinearity check
  svd.res <- tolerance_svd(ours_scale(data,center=center,scale=scale), tol = NULL)
  if( any(svd.res$d^2 < tol) ) {
    if(collinearity.stop){
      stop("continuous_mcd: Matrix is singular. It is likely that data are collinear where some variables are combinations of other variables.")
    }
  }
  ## alpha check
  if(alpha < .5){
    alpha <- .5
  }

  ## sample finder
  mcd.samples <- continuous_mcd_find_sample(data, center=center, scale=scale, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent,tol=tol)

  ## get sample distances
  tsvd.res <- tolerance_svd(ours_scale(data,center=center,scale=scale),tol = tol)
  mahals <- rowSums(tsvd.res$u^2)
  chis <- rowSums( (tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T))^2 )

  ## get robust mean & cov (loadings)
  rob.sample <- ours_scale(data[mcd.samples$final.orders[1,],],center=center,scale=scale)
  rob.center <- attributes(rob.sample)$`scaled:center`
  rob.scale <- attributes(rob.sample)$`scaled:scale`
  robust.tsvd.res <- tolerance_svd(rob.sample,tol = tol)

  ## get robust distances
  robust.vectors.and.scores <- continuous_scores_dists(data,center=rob.center,scale=rob.scale,robust.tsvd.res$v,robust.tsvd.res$d)
  # robust.mahals <- rowSums(robust.vectors.and.scores$sup.u^2)
  # robust.chis <- rowSums(robust.vectors.and.scores$sup.fi^2)


  res <- list(
    cov = list(loadings = robust.tsvd.res$v,
               singular.values = robust.tsvd.res$d,
               center = rob.center,
               scale = rob.scale
    ),
    dists = list(rob.md = robust.vectors.and.scores$mahals,
                 rob.chid = robust.vectors.and.scores$chis,
                 md = mahals,
                 chid = chis),
    det.samps = list(dets = mcd.samples$final.dets,
                     samples = mcd.samples$final.orders)
  )
  class(res) <- c("list", "contMCD")


  res
}









#'
#'  @export

continuous_mcd_find_sample <- function(data, center=T, scale=F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  ## data size check
  if(ncol(data) > (nrow(data)*.9)){
    stop("continuous_mcd: the column-to-row ratio is too high so 'mcd' cannot be performed")
  }
  ## alpha check
  if(alpha < .5){
    alpha <- .5
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
      init.norm <- ours_scale(data[init.samp,],center,scale)
      init.svd <- tolerance_svd(init.norm, tol = tol)
      init.mds <- round(rowSums(init.svd$u^2),digits=8)

      if(length(unique(init.mds)) < 2){
        init.size <- init.size + 1
      }else{
        sup.scores.and.dists <- continuous_scores_dists(data, attributes(init.norm)$`scaled:center`, attributes(init.norm)$`scaled:scale`, init.svd$v, init.svd$d)
        samp.config <- sort(order(sup.scores.and.dists$mahals)[1:h.size])
        findInit <- F
      }
    }

    min.info <- continuous_c_step(data, samp.config, center, scale, max.det.iters)
    dets[i] <- min.info$min.det
    orders[i,] <- min.info$obs.order
  }

  perc.cut <- round(num.subsets * top.sets.percent)
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique(unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),])

  final.dets <- vector("numeric", nrow(final.configs))
  final.orders <- matrix(NA,nrow(final.configs),h.size)
  for( i in 1:nrow(final.configs)){

    min.info <- continuous_c_step(data, final.configs[i,], center, scale, 100)	## set to Inf so that this converges on its own; need to make this settable & have a real max embedded in c.step
    final.dets[i] <- min.info$min.det
    final.orders[i,] <- min.info$obs.order

  }

  best.order <- order(final.dets)
  final.dets <- final.dets[best.order]
  final.orders <- final.orders[best.order,]

  return( list(final.dets = final.dets, final.orders = final.orders) )

}



#'
#' @export

continuous_scores_dists <- function(data,center=F,scale=F,loadings,singular.values){

  sup.fi <- (ours_scale(data,center=center,scale=scale) %*% loadings)
  sup.u <- sweep(sup.fi,2,singular.values,"/")

  mahals <- rowSums(sup.u^2)
  chis <- rowSums(sup.fi^2)

  return( list(mahals = mahals, chis = chis, sup.u=sup.u, sup.fi=sup.fi) )

}





#'
#' @export
#'


continuous_c_step <- function(data, obs.order, center=T, scale=F, max.iters=25, tol=sqrt(.Machine$double.eps)){

  old.det <- Inf
  old.center <- NaN
  old.v <- matrix(NaN,0,0)
  new.order <- old.order <- obs.order

  for(i in 1:max.iters){

    sub.data <- data[new.order,]
    sub.data.normed <- ours_scale(sub.data,center=center,scale=scale)
    svd.res <- tolerance_svd(sub.data.normed, tol = tol)

    new.center <- attributes(sub.data.normed)$`scaled:center`
    new.scale <- attributes(sub.data.normed)$`scaled:scale`
    new.det <- geometric.mean(svd.res$d^2)

    if( (new.det <= old.det) & (!isTRUE(all.equal(new.det,0,tolerance=tol))) ){
      if( center.sigma_checker(old.center, new.center, old.v, svd.res$v,tol=tol) & isTRUE(all.equal(new.det, old.det, tolerance= tol)) ){
        return( list(obs.order = old.order, min.det = old.det) )
      }
      else{

        old.det <- new.det
        old.center <- new.center
        old.v <- svd.res$v
        old.order <- new.order

        sup.scores.and.dists <- continuous_scores_dists(data, new.center, new.scale, svd.res$v, svd.res$d)
        new.order <- sort(order(sup.scores.and.dists$mahals)[1:length(obs.order)])

        if( isTRUE(all.equal(sort(new.order),sort(old.order))) ){
          return( list(obs.order = old.order, min.det = old.det) )
        }
      }
    }else{
      return( list(obs.order = old.order, min.det = old.det) )
    }
  }

  return( list(obs.order = new.order, min.det= new.det) )
}

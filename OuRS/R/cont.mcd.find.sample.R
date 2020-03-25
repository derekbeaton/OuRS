
## I think maybe I do not want people to interface directly with this function.
## that way I don't need to perform the checks here, just the checks in continuous_mcd
## or I can make this available but I'll need tests here, too...
### minimal tests here.. but those only override stupid things. also if you use this, you're committed to being OK with singular matrices
#'
#'  @export

cont_mcd_find_sample <- function(data, center=T, scale=F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  ## rudimentary data checks, or no?

  ## check alpha and change it if it's a problem.
  if(alpha<.5){
    # h.size <- floor((nrow(data)+1)/2)
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
        sup.scores <- cont.sup.fi.u(data, attributes(init.norm)$`scaled:center`, attributes(init.norm)$`scaled:scale`, init.svd$v, init.svd$d)
        mahals <- rowSums(sup.scores$sup.u^2)
        samp.config <- sort(order(mahals)[1:h.size])
        findInit <- F
      }
    }

    min.info <- cont.c.step(data, samp.config, center, scale, max.det.iters)
    dets[i] <- min.info$min.det
    orders[i,] <- min.info$obs.order
  }

  perc.cut <- round(num.subsets * top.sets.percent)
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique(unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),])

  final.dets <- vector("numeric", nrow(final.configs))
  final.orders <- matrix(NA,nrow(final.configs),h.size)
  for( i in 1:nrow(final.configs)){

    min.info <- cont.c.step(data, final.configs[i,], center, scale, 100)	## set to Inf so that this converges on its own; need to make this settable & have a real max embedded in c.step
    final.dets[i] <- min.info$min.det
    final.orders[i,] <- min.info$obs.order

  }

  best.order <- order(final.dets)
  final.dets <- final.dets[best.order]
  final.orders <- final.orders[best.order,]

  return( list(final.dets = final.dets, final.orders = final.orders) )

}

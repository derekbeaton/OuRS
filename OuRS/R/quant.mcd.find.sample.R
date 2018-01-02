

### the point of this function is to strictly return the best sample(s) for MCD.
  ## however, this is the heavy-duty part of MCD.

quant.mcd.find.sample <- function(data, center=T, scale=F, alpha=.75, h.size.abs=F, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05){


  if(h.size.abs){
    if(alpha<=.5){
      h.size <- floor((nrow(data)+1)/2)
    }else{
      h.size <- ceiling(nrow(data)*alpha)
    }
  }else{
    h.size <- h.alpha.n(alpha,nrow(data),ncol(data))
  }
  max.det.iters <- round(max.total.iters / num.subsets)

  dets <- vector("numeric", num.subsets)
  orders <- matrix(NA,num.subsets,h.size)
  for(i in 1:num.subsets){

    ##actually, make a while here to find the best random min(dim(data)+1 set, then pass the order in
    ##instead of det = 0, I can test if the mahals are all ~= from the min(dim(data)+1 set
    findInit <- T
    ## this doesn't make much sense... this needs to be smarter.
    init.size <- min(dim(data))+1

    while( findInit ){

      init.samp <- sort(sample(nrow(data),init.size))
      init.norm <- expo.scale(data[init.samp,],center,scale)
      init.svd <- tolerance.svd(init.norm)
      init.mds <- round(rowSums(init.svd$u^2),digits=8)

      if(length(unique(init.mds)) < 2){
        init.size <- init.size + 1
      }else{
        sup.scores <- sup.fi.u(data, attributes(init.norm)$`scaled:center`, attributes(init.norm)$`scaled:scale`, init.svd$v, init.svd$d)
        mahals <- rowSums(sup.scores$sup.u^2)
        samp.config <- sort(order(mahals)[1:h.size])
        findInit <- F
      }
    }

    min.info <- quant.c.step(data, samp.config, center, scale, max.det.iters)	## set to some threshold to not go for a long time
    dets[i] <- min.info$min.det
    orders[i,] <- min.info$obs.order
  }

  perc.cut <- round(num.subsets * top.sets.percent)	## they take "top 10" -- I will allow our search to be broader
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),]


  final.dets <- vector("numeric", nrow(final.configs))
  final.orders <- matrix(NA,nrow(final.configs),h.size)
  for( i in 1:nrow(final.configs)){

    ## RETURN TO THIS. I NEED INF to work.
    min.info <- quant.c.step(data, final.configs[i,], center, scale, 1000)	## set to Inf so that this converges on its own; need to make this settable & have a real max embedded in c.step
    final.dets[i] <- min.info$min.det
    final.orders[i,] <- min.info$obs.order
  }
  best.order <- order(final.dets)
  final.dets <- final.dets[best.order]
  final.orders <- final.orders[best.order,]


  return( list(final.dets = final.dets, final.orders = final.orders) )

}



### the point of this function is to strictly return the best sample(s) for MCD.
  ## however, this is the heavy-duty part of MCD.

cat.mcd.find.sample <- function(data, make.data.disjunctive=F,alpha=.75,num.subsets=500,max.total.iters=num.subsets*20,top.sets.percent=.05){
  if(make.data.disjunctive){
    data <- make.data.nominal(data)
  }
  preproc.data <- ca.preproc(data) ## this could be more efficient...
  profiles <- (diag(1/preproc.data$m) %*% preproc.data$Ox)

  if(alpha<.5){
    h.size <- floor((nrow(data)+1)/2)
  }
  h.size <- h.alpha.n(alpha,nrow(data),ncol(data))
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

      init.svd <- tolerance.svd(preproc.data$weightedZx[init.samp,])
      init.mds <- round(rowSums(init.svd$u^2),digits=8)	## do I need to round?

      if(length(unique(init.mds)) < 2){
        init.size <- init.size + 1
      }else{
        sup.scores <- mahal.from.ca(profiles, preproc.data$m, preproc.data$w, init.svd$v, init.svd$d) ## this really needs updating.

        samp.config <- sort(order(sup.scores)[1:h.size])
        findInit <- F
      }
    }

    min.info <- cat.c.step(profiles, preproc.data$weightedZx, preproc.data$m, preproc.data$w, samp.config, max.det.iters)
    dets[i] <- min.info$min.det
    orders[i,] <- min.info$obs.order
  }

  perc.cut <- round(num.subsets * top.sets.percent)	## they take "top 10" -- I will allow our search to be broader
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),]
  ### I should add in a best and worst case order for final.configs


  final.configs <- unique(final.configs)


  final.dets <- vector("numeric", nrow(final.configs))
  final.orders <- matrix(NA,nrow(final.configs),h.size)
  for( i in 1:nrow(final.configs)){
    min.info <- cat.c.step(profiles, preproc.data$weightedZx, preproc.data$m, preproc.data$w, eigen.fix=eigen.fix, final.configs[i,], 1000)		## set to Inf so that this converges on its own; need to make this settable & have a real max embedded in c.step
    final.dets[i] <- min.info$min.det
    final.orders[i,] <- min.info$obs.order
  }

  best.order <- order(final.dets)
  final.dets <- final.dets[best.order]
  final.orders <- final.orders[best.order,]
  return( list(final.dets= final.dets, final.orders= final.orders) )
}

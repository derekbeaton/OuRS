cat.mcd <- function(data, make.data.disjunctive=F, alpha=.75, h.size.abs=F, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  if(make.data.disjunctive){
    data <- make.data.nominal(data)
  }

  if(ncol(data) > (nrow(data)*.9)){
    stop("Column:Row ratio is too high. Please use another method (e.g., robPCA, rPCA, POET).")
  }


  ## sample finder
  mcd.samples <- cat.mcd.find.sample(data, make.data.disjunctive=make.data.disjunctive, alpha=alpha, h.size.abs=h.size.abs, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent)
  ## only grab the top sample.
  best.sample <- mcd.samples$final.orders[1,]


  ### NOTE: I should pass the tol parameter through to all places where tolerance.svd gets called.
  ## compute standard distances (especially Mahal)
  # tsvd.res <- tolerance.svd(expo.scale(data,center=center,scale=scale))
  # mahals <- rowSums(tsvd.res$u^2)
  # chis <- rowSums( (tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T))^2 )



}

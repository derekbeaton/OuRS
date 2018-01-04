cat.c.step <- function(profiles, weighted.deviations, row.weights, col.weights, obs.order, max.iters=25, tol=sqrt(.Machine$double.eps)){

  old.det <- Inf
  old.center <- NaN
  old.v <- matrix(NaN,0,0)
  new.order <- old.order <- obs.order

  for(i in 1:max.iters){
    sub.data <- weighted.deviations[new.order,]
    new.center <- colMeans(sub.data)
    svd.res <- tryCatch( {tolerance.svd(sub.data)}, error=function(x) 'FAIL') ## this should probably be explained a bit!

    if(length(svd.res)==3){
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

          sup.scores <- cat.sup.fi.u(profiles, row.weights, col.weights, svd.res$v, svd.res$d)
          mahals <- rowSums(sup.scores$sup.u^2)
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

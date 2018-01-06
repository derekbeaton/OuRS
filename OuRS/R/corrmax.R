

#cont.corrmax <- function(population.data,sample.data){
cont.corrmax <- function(target.data,rob.center=T,rob.scale=F,loadings,singular.values,tol=.Machine$double.eps){

  ## if we want to do a correction for sample size in this case, all we need to do is scale up the sv/eigen values.

  diag.sampcov <- sqrt(diag(tcrossprod( sweep(loadings,2,singular.values,"*") )))
  inv.DSD.half <- (tcrossprod(sweep( sweep(loadings,2,singular.values,"*"),1,diag.sampcov,"/") %^% (-1/2)))
  target.data <- expo.scale(target.data,center=rob.center,scale=rob.scale) #signs are switched; I can just *-1
  W <- target.data %*% sweep(inv.DSD.half,1,diag.sampcov,"/")
  w.svd <- tolerance.svd(W,tol=tol)

    ## we should use percentages in place of contributions as they have meaning
    ## also: we now have a new Mahalanobis from corrmax.
  return(list(mah=rowSums(w.svd$u^2),percs=sweep(W^2,1,rowSums(W^2),"/")*100))

}



## this needs to be expanded / updated
## I need a plain corrmax and a reference corrmax (which is what this is)
## and I need to return more than just W^2


#cont.corrmax <- function(population.data,sample.data){
continous_corrmax <- function(target.data,rob.center=T,rob.scale=F,loadings,singular.values,tol=.Machine$double.eps){


  target.data <- ours_scale(target.data,center=rob.center,scale=rob.scale)

  diag.sampcov <- sqrt(diag(tcrossprod( sweep(loadings,2,singular.values,"*") )))

  inv.DSD.half <- tcrossprod( GSVD::invsqrt_psd_matrix(sweep( sweep(loadings,2,singular.values,"*"),1,diag.sampcov,"/")) )

  W <- target.data %*% sweep(inv.DSD.half,1,diag.sampcov,"/")
  w.svd <- tolerance_svd(W,tol=tol)
    ## no need because, as in Eq 2(?) in Garthwaite & Koch (2016) they point out that X --> W exists where mahal(X) == mahal(W)
    ## but I've brought it back for convenience

    ## we should use percentages in place of contributions as they have meaning
  #return(list(mah=rowSums(w.svd$u^2),percs=sweep(W^2,1,rowSums(W^2),"/")*100))
  
  res <- list(
    mahal_dists = rowSums(w.svd$u^2),
    mahal_dists2 = diag(crossprod(W)),
    corr_max_transform = inv.DSD.half,
    transformed_data = W,
    percentage_contributions = sweep(W^2,1,rowSums(W^2),"/")*100
  )
  
  class(res) <- c("list", "OuRS", "CorrMax", "continuous")
  return(res)

}


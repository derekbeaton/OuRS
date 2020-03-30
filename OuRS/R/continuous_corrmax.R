
## this needs to be expanded / updated
## I need a plain corrmax and a reference corrmax (which is what this is)
## and I need to return more than just W^2


#cont.corrmax <- function(population.data,sample.data){

#' @title CorrMax transformation
#' 
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#' 
#' @param target.data
#' @param rob.center
#' 
continous_corrmax <- function(target.data, center=T, scale=F, loadings, singular.values, tol=.Machine$double.eps){

  target.data <- ours_scale(target.data, center=center,scale=scale)
  
  if(missing(loadings) | missing(singular.values)){
    
    svd.res <- tolerance_svd(target.data, tol=tol)
    loadings <- svd.res$v
    singular.values <- svd.res$d
    svd.res <- NULL ## get rid of it
    
  }
  
  diag.sampcov <- sqrt(diag(tcrossprod( sweep(loadings,2,singular.values,"*") )))
    ## this can be made more efficient if I use /svd()$d as another sweep inside the tcrossprod()
  inv.DSD.half <- GSVD::invsqrt_psd_matrix( tcrossprod( sweep( sweep(loadings,2,singular.values,"*"),1,diag.sampcov,"/") ) )
  W <- target.data %*% sweep(inv.DSD.half,1,diag.sampcov,"/")
  
  res <- list(
    corr_max_transform = inv.DSD.half,
    transformed_data = W,
    mahal_dists = rowSums(W^2),
    percentage_contributions = sweep(W^2,1,rowSums(W^2),"/")*100
  )
  
  class(res) <- c("list", "OuRS", "CorrMax", "continuous")
  return(res)

}


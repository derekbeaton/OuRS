
## this needs to be expanded / updated
## I need a plain corrmax and a reference corrmax (which is what this is)
## and I need to return more than just W^2


## I need to rename some things ("target.data") and also work out the needs of the generalized_ and friends approaches
  #### I need to see existing pipeline code because those may use profiles or weightedZx... I don't recall right now
  #### previous note indicates "## here target.data has to be the weighted deviations (weightedZx)" which makes sense.
    ### I suppose it *could* be Zx too, but the weightedZx is where we get the MD anyways.

#cont.corrmax <- function(population.data,sample.data){

#' @title CorrMax transformation
#' 
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#' 
#' @param target.data
#' @param center
#' @param scale
#' @param loadings
#' @param singular.values
#' 
#' @example 
#' 
#' @author Derek Beaton
#' @export
#' 
continous_corrmax <- function(target.data, center=T, scale=F, loadings, singular.values){

  target.data <- ours_scale(target.data, center=center,scale=scale)
  
  if(missing(loadings) | missing(singular.values)){
    
    pca.res <- pca(target.data, tol=tol)
    loadings <- pca.res$v
    singular.values <- pca.res$d
    pca.res <- NULL ## get rid of it
    
  }
  
  res <- corrmax_core_tranform(target.data, loadings, singular.values)
  
  class(res) <- append(class(res), "continuous")
  return(res)

}


categorical_corrmax <- function(target.data,loadings,singular.values,tol=.Machine$double.eps){
  
  ## if we want to do a correction for sample size in this case, all we need to do is scale up the sv/eigen values.
  
  #target.data <- expo.scale(target.data,center=rob.center,scale=rob.scale) #signs are switched; I can just *-1
  
  # diag.sampcov <- sqrt(diag(tcrossprod( sweep(loadings,2,singular.values,"*") )))
  # inv.DSD.half <- (tcrossprod(sweep( sweep(loadings,2,singular.values,"*"),1,diag.sampcov,"/") %^% (-1/2)))
  # 
  # W <- target.data %*% sweep(inv.DSD.half,1,diag.sampcov,"/")
  # 
  # return(sweep(W^2,1,rowSums(W^2),"/")*100)

  ### need all that stuff here...
  
  ### actually, this should call a "generalized" corrmax I think...
  ### and this function here, like with the MCD, just passes stuff along after it does disjunctive_coding
  # res <- corrmax_core_tranform(target.data, loadings, singular.values)
  # 
  # class(res) <- append(class(res), "continuous")
  # return(res)
  
}


# ordinal_corrmx <- function(){
#   
# }

# generalized_corrmax <- function(){
#
# }


# mixed_data_corrmax <- function(){
#   
# }


#' @noRd
corrmax_core_tranform <- function(target.data, loadings, singular.values){
  
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
  
  class(res) <- c("list", "OuRS", "CorrMax")
  return(res)
  
}
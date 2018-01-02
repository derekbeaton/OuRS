svd.collinearity.test <- function(singular.values,tol=.Machine$double.eps){

  d.below.tol <- singular.values < tol
  collinear.components <- which(d.below.tol)
  any.d.below.tol <- any(d.below.tol)
  if( any.d.below.tol ){
    warn(
      paste0("Some singular values below tolerance. Several measures likely collinear. Please check components: ",collinear.components)
    )
  }
  return(collinear.components) #should be length > 0 if TRUE, length = 0 if F
}

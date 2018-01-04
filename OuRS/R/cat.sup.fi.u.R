cat.sup.fi.u <- function(profiles, row.weights, col.weights, loadings, singular.values){

  ## NOT EFFICIENT. MAKE MORE EFFICIENT

  sup.fi <- (profiles %*% (diag(sqrt(col.weights)*(1/col.weights)) %*% loadings))
  sup.u <- diag(sqrt(row.weights)) %*% (sup.fi * matrix(1/singular.values,nrow(profiles),ncol(loadings),byrow=T))
  return( list( sup.fi=sup.fi,sup.u=sup.u ) )

}

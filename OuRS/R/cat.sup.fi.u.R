cat.sup.fi.u <- function(profiles, row.weights, col.weights, loadings, singular.values){

  ## NOT EFFICIENT. MAKE MORE EFFICIENT
  ## also this can handle the profiles here...

  #sup.fi <- (profiles %*% (diag(sqrt(col.weights)*(1/col.weights)) %*% loadings))
  sup.fi <- profiles %*% sweep(loadings,1,sqrt(col.weights)/col.weights,"*")

  #sup.u <- diag(sqrt(row.weights)) %*% (sup.fi * matrix(1/singular.values,nrow(profiles),ncol(loadings),byrow=T))
  sup.u <- sweep(sweep(sup.fi,2,singular.values,"/"),1,sqrt(row.weights),"*")

  return( list( sup.fi=sup.fi,sup.u=sup.u ) )

}

## should this take in a class to determine which route to go?
  ## maybe for now (easy) just make these

cont.sup.fi.u <- function(data,center=F,scale=F,loadings,singular.values){

  sup.fi <- (expo.scale(data,center=center,scale=scale) %*% loadings)
  sup.u <- sup.fi * matrix(1/singular.values,nrow(data),ncol(loadings),byrow=T)

  return( list(sup.fi=sup.fi,sup.u=sup.u) )

}

## should this take in a class to determine which route to go?
  ## maybe for now (easy) just make these

sup.fi.u <- function(data,center=F,scale=F,loadings,singular.values){

  fi.scores <- (expo.scale(data,center=center,scale=scale) %*% loadings)
  u.scores <- fi.scores * matrix(1/singular.values,nrow(data),ncol(loadings),byrow=T)

  return( list(sup.fi=fi.scores,sup.u=u.scores) )

}

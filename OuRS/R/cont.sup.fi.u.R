## this should also directly return scores.

cont.sup.fi.u <- function(data,center=F,scale=F,loadings,singular.values){

  sup.fi <- (ours_scale(data,center=center,scale=scale) %*% loadings)
  sup.u <- sweep(sup.fi,2,singular.values,"/")
  return( list(sup.fi=sup.fi,sup.u=sup.u) )

}


cont.boot.sup.fi.u <- function(target.data,center=T,scale=F,loadings,singular.values,iters=100){

  boot.distrs <- matrix(NA,nrow(target.data),iters)
  for(i in 1:iters){

    boot.distrs[,i] <- rowSums( (expo.scale(target.data[sample(nrow(target.data),replace=T),],center=center,scale=scale) %*% loadings %*% diag(1/singular.values))^2 )

  }
  return(boot.distrs)
}

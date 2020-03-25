
cont.boot.sup.u <- function(target.data,center=T,scale=F,loadings,singular.values,iters=100){

  boot.distrs <- matrix(NA,nrow(target.data),iters)
  for(i in 1:iters){

    boot.distrs[,i] <- rowSums(

      sweep( (ours_scale(target.data[sample(nrow(target.data),replace=T),],center=center,scale=scale) %*% loadings) , 2, singular.values, "/")^2

      )

  }
  return(boot.distrs)
}

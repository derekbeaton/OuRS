
cat.boot.sup.fi.u <- function(target.data,loadings,singular.values,iters=100){

  boot.distrs <- matrix(NA,nrow(target.data),iters)
  for(i in 1:iters){
    ca.preproc.data <- ca.preproc(target.data[sample(nrow(target.data),replace=T),])

      ## I'm sure this could be much more efficient in some way.
    these.vecs <- sweep(loadings,1,sqrt(ca.preproc.data$w)/ca.preproc.data$w,"*")
    these.vecs[is.nan(these.vecs)] <- 0
    boot.distrs[,i] <- rowSums(
      (sweep(
        sweep(

            sweep(ca.preproc.data$Ox,1,ca.preproc.data$m,"/") %*% these.vecs

        ,2,singular.values,"/")
      ,1,sqrt(ca.preproc.data$m),"*"))^2
    )
  }
  return(boot.distrs)
}





ca <- function(X,k=0){

  sum.data <- sum(X)
  rowSums.data <- rowSums(X)
  wi <- rowSums.data/sum.data
  wj <- colSums(X)/sum.data
  ca.res <- gsvd( sweep(sweep(X,1,rowSums.data,"/"),2,wj), wi, 1/wj, k = k )
  ca.res$fi <- sweep(ca.res$fi,1,wi,"/")

  return( ca.res )
}

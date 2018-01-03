ca.preproc <- function(X){
  Ox <- X/sum(X)
  m <- rowSums(Ox)
  w <- colSums(Ox)
  Ex <- m %o% w
  Zx <- Ox - Ex
  weightedZx <- (sqrt.mat(diag(1/m)) %*% Zx %*% sqrt.mat(diag(1/w)))
  return( list(m=m,w=w,Zx=Zx,Ox=Ox,Ex=Ex,weightedZx=weightedZx) )
}

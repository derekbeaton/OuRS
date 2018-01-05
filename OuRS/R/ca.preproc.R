ca.preproc <- function(X){
  Ox <- X/sum(X)
  m <- rowSums(Ox)
  w <- colSums(Ox)
  Ex <- m %o% w
  Zx <- Ox - Ex
  #weightedZx <- (sqrt.mat(diag(1/m)) %*% Zx %*% sqrt.mat(diag(1/w)))
  #weightedZx <- ((diag(1/m) %^% (1/2)) %*% Zx %*% (diag(1/w) %^% (1/2)))
  weightedZx <- sweep(sweep(Zx,1,sqrt(m),"/"),2,sqrt(w),"/")
  return( list(m=m,w=w,Zx=Zx,Ox=Ox,Ex=Ex,weightedZx=weightedZx) )
}

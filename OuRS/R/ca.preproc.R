  ## obviously this causes a major memory hold...
    ## OX, EX, ZX, and weightedZx are all the same size as X...
    ## which means there are at least 5 matrices of identical size at any given time.
    ## I need to go back through the code to identify places where we really don't need all this.

ca.preproc <- function(X){
  Ox <- X/sum(X)
  m <- rowSums(Ox)
  w <- colSums(Ox)
  Ex <- m %o% w
  Zx <- Ox - Ex
  weightedZx <- sweep(sweep(Zx,1,sqrt(m),"/"),2,sqrt(w),"/")
  return( list(m=m,w=w,Zx=Zx,Ox=Ox,Ex=Ex,weightedZx=weightedZx) )
}

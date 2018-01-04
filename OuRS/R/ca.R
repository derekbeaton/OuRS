## test change.

ca <- function(X,k=0){
  preproc <- ca.preproc(X)
  ca.res <- gsvd(preproc$Zx, 1/preproc$m, 1/preproc$w, k=k)
  #ca.res$fi <- M %*% ca.res$p %*% diag(ca.res$Dv)
  #rownames(ca.res$fi) <- rownames(preproc$Zx)
  #ca.res$fj <- W %*% ca.res$q %*% diag(ca.res$Dv)
  #rownames(ca.res$fj) <- colnames(preproc$Zx)
  return( ca.res )
}

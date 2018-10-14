## Candes's.
  ## Effectively rewritten as a combination of: https://github.com/dlaptev/RobustPCA/blob/master/RobustPCA.m & 'rpca' package


my.rpca <- function(DATA, lambda = 1 / sqrt(max(dim(DATA))), mu = 10 * lambda, maximum.iterations = 1000, tol=sqrt(.Machine$double.eps)){

  shrinkage.operator <- function(X, tau){

    return( sign(X) * pmax( abs(X) - tau, 0 ) )

  }

  stopping.condition <- tol * sqrt(sum(DATA^2))

## intialize
  DATA <- as.matrix(DATA)
  L <- S <- Y <- matrix(0, nrow(DATA), ncol(DATA))

  for(i in 1:maximum.iterations){

    res <- tolerance.svd(DATA - S + Y, tol=tol)
    L <- sweep(res$u, 2, shrinkage.operator(res$d, 1/mu), "*") %*% t(res$v)

    S <- shrinkage.operator( DATA - L + Y, lambda/mu)
    Z <- DATA - L - S

    if(sqrt(sum(Z^2)) < stopping.condition){
      break
    }
    Y <- Y + Z
  }
  return(list(L=L, S=S, i=i))
}

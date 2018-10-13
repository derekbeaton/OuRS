## Candes's.
  ## Effectively rewritten as a combination of: https://github.com/dlaptev/RobustPCA/blob/master/RobustPCA.m & 'rpca' package


rpca <- function(DATA, lambda = 1 / sqrt(max(dim(DATA))), mu = 10 * lambda, maximum.iterations = 1000, tol=sqrt(.Machine$double.eps)){#, switcher = T){

  shrinkage.operator <- function(X, tau){

    return( sign(X) * pmax( abs(X) - tau, 0 ) )

  }
  reg.svd <- function(X, tau, tol=tol=sqrt(.Machine$double.eps)){
    res <- tolerance.svd(X, tol=tol)
    return( sweep(res$u, 2, shrinkage.operator(res$d, tau), "*") %*% t(res$v) )
  }

  stopping.condition <- tol * sqrt(sum(DATA^2))

## intialize
  DATA <- as.matrix(DATA)
  L <- S <- Y <- matrix(0, nrow(DATA), ncol(DATA))

  for(i in 1:maximum.iterations){
    # if(!switcher){
    #   L <- reg.svd(DATA - S + ((1/mu)*Y), 1/mu)
    # }else{
      L <- reg.svd(DATA - S + Y, 1/mu)
    #}

    S <- shrinkage.operator( DATA - L + Y, lambda/mu)
    Z <- DATA - L - S

    if(sum(sqrt(Z)) < stopping.condition){
      break
    }

    # if(!switcher){
    #   Y <- Y + mu*Z
    # }else{
      Y <- Y + Z
    #}
  }

}

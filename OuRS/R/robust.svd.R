## Candes's.
  ## Effectively rewritten as a combination of: https://github.com/dlaptev/RobustPCA/blob/master/RobustPCA.m & 'rpca' package


robust.svd <- function(DATA, lambda = 1 / sqrt(max(dim(DATA))), mu = 10 * lambda, maximum.iterations = 1000, tol=sqrt(.Machine$double.eps)){

  shrinkage.operator <- function(X, tau){

    return( sign(X) * pmax( abs(X) - tau, 0 ) )

    ## from Verbanck's regularized PCA; but we need to know dimensionality. Perhaps for later. Or it can go elsewhere.
    # sigma2 = sum(res$d[-c(1:S)]^2)/(n*p - p - n*S - p*S + S^2)
    # lambda.shrinked = (res$d[1:S]^2 - n*(p/min(p, (n-1)))*sigma2)/res$d[1:S]

  }


  DATA <- as.matrix(DATA)

  stopping.condition <- tol * sqrt(sum(DATA^2))
  L <- S <- Y <- matrix(0, nrow(DATA), ncol(DATA))

  for(i in 1:maximum.iterations){

    res <- tolerance.svd(DATA - S + Y, tol=tol)
    L <- sweep(res$u, 2, shrinkage.operator(res$d, 1/mu), "*") %*% t(res$v)
    S <- shrinkage.operator( DATA - L + Y, lambda/mu)


    if(sqrt(sum( (DATA - L - S) ^2)) < stopping.condition){
      break
    }
    Y <- Y + (DATA - L - S)
  }
  return(list(L=L, S=S, i=i))
}

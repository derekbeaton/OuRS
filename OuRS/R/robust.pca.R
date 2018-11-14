robust.pca <- function(DATA, center=T, scale=T, lambda = 1 / sqrt(max(dim(DATA))), mu = 10 * lambda, maximum.iterations = 1000, tol=sqrt(.Machine$double.eps)){

  DATA <- expo.scale(DATA)
    ## attributes are retained.

}

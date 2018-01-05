## Modified from: https://stackoverflow.com/questions/11057639/identifying-if-only-the-diagonal-is-non-zero
is.diagonal.matrix <- function(x,tol=.Machine$double.eps){
  if(is.null(dim(x))){
    stop("is.diagonal.matrix: X is not a matrix.")
  }
  x[ x^2 < tol ] <- 0
  return(all(x[lower.tri(x)] == 0, x[upper.tri(x)] == 0))
}

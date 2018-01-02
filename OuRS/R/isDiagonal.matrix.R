#'
#'  @export
#'
#'  @title \code{isDiagonal.matrix}: test if a matrix is diagonal
#'
#'  @description \code{isDiagonal.matrix} takes a matrix and tests if it is diagonal matrix.
#'
#'  @param X a matrix to test
#'
#'  @return
#'  boolean. TRUE if the matrix is diagonal, FALSE if the matrix is not.


  ## I stole this from somewhere... but I don't remember where
isDiagonal.matrix <- function(X){
  if(is.null(dim(X))){
    stop("isDiagonal.matrix: X is not a matrix.")
  }
  return(all(X[lower.tri(X)] == 0, X[upper.tri(X)] == 0))
}

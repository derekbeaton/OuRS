## stolen from rrcov or robustbase
#' @title H subset from alpha
#' @description Compute the size of the H subset from alpha, number of rows (n), and number of columns (p)
#' 
#' @param alpha numeric. Percentage value for computing subsample size. Should exist between [.5,1]. If outside of that range it will be changed to .5 or 1.
#' @param n number of rows for data to compute the H subset
#' @param p number of columns for data to compute the H subset
#' @noRd

h.alpha.n <- function (alpha, n, p)
{
  if(alpha < .5){
    alpha <- .5
  }
  if(alpha > 1){
    alpha <- 1
  }
  n2 <- (n + p + 1)%/%2
  floor(2 * n2 - n + 2 * (n - n2) * alpha)
}




### this is a check to see if the new center and covariance are identical to the old, so we can stop the search
#' @title Check old and new data center and loadings matrices
#' @description a check that returns a logical. If the center and loadings matrix of iteration \eqn{(n-1)} is the same as center and loadings matrix of iteration \eqn{n}, then return \code{TRUE}. Else \code{FALSE}
#' @details This function is meant to test if the MCD search should stop because the current iteration has the same center and loadings as the previous iteration.
#' @param old.center a numeric vector. The column-wise center from iteration \eqn{n-1}
#' @param new.center a numeric vector. The column-wise center from the current iteration \eqn{n}
#' @param old.v a numeric matrix. The column loadings (from the SVD) of a matrix from iteration \eqn{n-1}
#' @param new.v a numeric matrix. The column loadings (from the SVD) of a matrix from the current iteration \eqn{n}
#' @param tol numeric scalar >= 0. Any values smaller than \code{tol} are considered 0
#' @noRd

center.sigma_checker <- function(old.center=NaN,new.center=NaN,old.v=matrix(NaN,0,0),new.v=matrix(NaN,0,0),tol=sqrt(.Machine$double.eps)){
  
  ## not entirely sure why these two if() statements are here.
  ### the center & span of columns should be the same.
  if(length(old.center)!=length(new.center)){
    old.center <- rep(0,length(new.center))
  }
  if( nrow(old.v)!=nrow(new.v) | ncol(old.v)!=ncol(new.v) ){
    old.v <- matrix(0,nrow(new.v),ncol(new.v))
  }
  ## just make sure that the differences between the two items are negligible
  if( (sum( sqrt((old.center - new.center)^2) ) < tol) & (sum( sqrt((old.v - new.v)^2) ) < tol) ){
    return(TRUE)
  }
  return(FALSE)
}



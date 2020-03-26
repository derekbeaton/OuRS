## don't export any of these, but do document them?


## stolen from 'psych' package
#' @export
#'

geometric_mean <- function (x, na.rm = TRUE)
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = TRUE))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}



## stolen from rrcov or robustbase
#' @export
#'

h.alpha.n <- function (alpha, n, p)
{
  n2 <- (n + p + 1)%/%2
  floor(2 * n2 - n + 2 * (n - n2) * alpha)
}




### this is a check to see if the new center and covariance are identical to the old, so we can stop the search

#' @export
#'
#'

center.sigma_checker <- function(old.center=NaN,new.center=NaN,old.v=matrix(NaN,0,0),new.v=matrix(NaN,0,0),tol=sqrt(.Machine$double.eps)){
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



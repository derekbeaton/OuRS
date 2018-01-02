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

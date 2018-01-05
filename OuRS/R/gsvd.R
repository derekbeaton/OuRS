#' @export
#'
#' @title Generalized SVD
#'
#' @description
#' \code{gsvd} takes in left (\code{LW}) and right (\code{RW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to the data (\code{DAT})
#'   Left and right constraints are used for the orthogonality conditions.
#'
#' @param DAT a data matrix to decompose
#' @param LW \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the matrix and thus left singular vectors.
#' @param RW \bold{R}ight \bold{W}eights -- the constraints applied to the right side (rows) of the matrix and thus right singular vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d.orig}). If 0, the full set of components are returned.
#' @param tol default is .Machine$double.eps. A parameter with two roles: A tolerance level for (1) eliminating (tiny variance or negative or imaginary) components and (2) converting all values < tol to 0 in \item{u} and \item{v}.
#'
#' @return A list with nine elements:
#' \item{d.orig}{A vector containing the singular values of DAT > \code{tol}.}
#' \item{tau}{A vector that contains the (original) explained variance per component (eigenvalues derived from \code{$d.orig}.}
#' \item{d}{A vector containing the singular values of x > \code{tol}. Length is \code{min(length(d.orig), k)}}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{p}{Left (rows) generalized singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{fi}{Left (rows) component scores. Dimensions are \code{nrow(DAT)} by k.}
#' \item{v}{Right (columns) singular vectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{q}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{fj}{Right (columns) component scores. Dimensions are \code{ncol(DAT)} by k.}
#'
#' @seealso \code{\link{tolerance.svd}} and \code{\link{svd}}
#'
#' @examples
#'  ## an example with correspondence analysis.
#'  data(authors)
#'  author.data <- authors$ca$data
#'  Observed <- author.data/sum(author.data)
#'  row.w <- rowSums(Observed)
#'    row.W <- diag(1/row.w)
#'  col.w <- colSums(Observed)
#'    col.W <- diag(1/col.w)
#'  Expected <- row.w %o% col.w
#'  Deviations <- Observed - Expected
#'  ca.res <- gsvd(Deviations,row.W,col.W)
#'
#'  ## an example with canonical correlation analysis (though not all bells and whistles exist)
#'  data(two.table.wine)
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'  cca.res <- gsvd(
#'      matrix.generalized.inverse(crossprod(X)) %*% t(X) %*% Y %*% matrix.generalized.inverse(crossprod(Y)),
#'      crossprod(X),
#'      crossprod(Y)
#'  )
#'
#'  cca.res$lx <- (X %*% cca.res$p)
#'  cca.res$ly <- (Y %*% cca.res$q)
#'
#'  \dontrun{
#'      optimize.for <- t(cca.res$lx) %*% cca.res$ly
#'      all.equal(diag(optimize.for),cca.res$d)
#'
#'
#'      base.cca <- cancor(X,Y,F,F)
#'
#'      sum(abs(base.cca$cor - cca.res$d)) < (.Machine$double.eps*100)
#'      base.cca$xcoef / cca.res$p
#'      base.cca$ycoef[,1:ncol(cca.res$q)] / cca.res$q
#'  }
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen


gsvd <- function(DAT, LW, RW, k = 0, tol=.Machine$double.eps){

  is.identity.matrix <- function(x,tol=.Machine$double.eps){
    if(is.null(dim(x))){
      stop("is.identity.matrix: x is not a matrix.")
    }
    x <- as.matrix(x)
    x[abs(x) < tol] <- 0
    if(is.diagonal.matrix(x)){
      if( all(diag(x)==1) ){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  }
  is.empty.matrix <- function(x,tol=.Machine$double.eps){
    x <- as.matrix(x)
    x[abs(x) < tol] <- 0
    if(sum(x)==0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }


  # preliminaries
  DAT.dims <- dim(DAT)
  if(length(DAT.dims)!=2){
    stop("gsvd: DAT must have dim length of 2 (i.e., rows and columns)")
  }
  DAT <- as.matrix(DAT)
  DAT[abs(DAT) < tol] <- 0
  RW.is.vector <- LW.is.vector <- RW.is.missing <- LW.is.missing <- F

  if(is.empty.matrix(LW)){
    stop("gsvd: LW is empty (i.e., all 0s")
  }
  if(is.empty.matrix(RW)){
    stop("gsvd: RW is empty (i.e., all 0s")
  }

  # check if LW and RW are missing, if they are vectors, or if they are diagonal matrices.

  if( missing(LW) ){
    LW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(LW) ) {
      LW.is.vector <- T
    }else if(!LW.is.vector){

      if( is.identity.matrix(LW) ){
        LW.is.missing <- T
        warning("gsvd: LW was an identity matrix. LW will not be used in the GSVD.")
      }else if( is.diagonal.matrix(LW) ){

        LW <- diag(LW)

        if( length(LW) != DAT.dims[1] ){
          stop("gsvd:length(LW) does not equal nrow(DAT)")
        }else{
          LW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(LW) != ncol(LW) | nrow(LW) != DAT.dims[1] ){
        stop("gsvd:nrow(LW) does not equal ncol(LW) or nrow(DAT)")
      }
    }
  }


  if( missing(RW) ){
    RW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(RW) ) {
      RW.is.vector <- T
    }else if(!RW.is.vector){

      if( is.identity.matrix(RW) ){
        RW.is.missing <- T
        warning("gsvd: RW was an identity matrix. RW will not be used in the GSVD.")
      }else if( is.diagonal.matrix(RW) ){

        RW <- diag(RW)

        if( length(RW) != DAT.dims[2] ){
          stop("gsvd:length(RW) does not equal ncol(DAT)")
        }else{
          RW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(RW) != ncol(RW) | nrow(RW) != DAT.dims[2] ){
        stop("gsvd:nrow(RW) does not equal ncol(RW) or ncol(DAT)")
      }
    }
  }


  ## these tests can be moved up but I just can't find a good place for them.
  if( LW.is.vector ){  ## replace with sweep
    DAT <- sweep(DAT,1,sqrt(LW),"*")
  }else if(!LW.is.missing){
    DAT <- (LW %^% (1/2)) %*% DAT
  }else{
    stop("gsvd: unknown condition for LW.")
  }

  if( RW.is.vector ){  ## replace with sweep
    DAT <- sweep(DAT,2,sqrt(RW),"*")
  }else if(!RW.is.missing){
    DAT <- DAT %*% (RW %^% (1/2))
  }else{
    stop("gsvd: unknown condition for RW.")
  }


  if(k<=0){
    k <- min(nrow(DAT),ncol(DAT))
  }

  res <- tolerance.svd(DAT,nu=k,nv=k,tol=tol)

  res$d.orig <- res$d
  res$tau <- res$d.orig^2/sum(res$d.orig^2)
  components.to.return <- min(length(res$d.orig),k) #a safety check

  res$d <- res$d.orig[1:components.to.return]
  ## u and v should already be k vectors but again, be safe.
  res$u <- as.matrix(res$u[,1:components.to.return])
  res$v <- as.matrix(res$v[,1:components.to.return])


  if(LW.is.vector){
    res$p <- sweep(res$u,1,1/sqrt(LW),"*")
    res$fi <- sweep(sweep(res$p,1,LW,"*"),2,res$d,"*")
  }else if(!LW.is.missing){
    res$p <- (LW %^% (-1/2)) %*% res$u
    res$fi <- sweep((LW %*% res$p),2,res$d,"*")
  }else{
    res$p <- res$u
    res$fi <- sweep(res$p,2,res$d,"*")
  }

  if(RW.is.vector){
    res$q <- sweep(res$v,1,1/sqrt(RW),"*")
    res$fj <- sweep(sweep(res$q,1,RW,"*"),2,res$d,"*")
  }else if(!RW.is.missing){
    res$q <- (RW %^% (-1/2)) %*% res$v
    res$fj <- sweep((RW %*% res$q),2,res$d,"*")
  }else{
    res$q <- res$v
    res$fj <- sweep(res$q,2,res$d,"*")
  }

  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- rownames(DAT)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(DAT)

  return(res)
  #return(list(fi = fi, fj = fj, p = p, q = q, u = res$u, v = res$v, d = d, d.orig = d.orig, tau = tau))
}

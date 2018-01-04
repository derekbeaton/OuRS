## a generalized SVD function.

#  Generalized SVD: A GSVD function that takes in left and right constraints (usually diagonal matrices, but any positive semi-definite matrix is fine).
#   Constraints are applied to the left and right singular vectors for the orthogonality constraint.

#'
#'  @export
#'
#'  @title \code{gsvd}: the generalized singular value decomposition.
#'
#'  @description \code{gsvd} takes in left (\code{LW}) and right (\code{RW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to the data (\code{DAT})
#'   Left and right constraints are used for the orthogonality conditions.
#'
#'  @param DAT a data matrix to decompose
#'  @param LW \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the matrix and thus left singular vectors
#'  @param RW \bold{R}ight \bold{W}eights -- the constraints applied to the right side (rows) of the matrix and thus right singular vectors
#'  @param nu the number of left singular vectors to be computed. Default is \code{min(dim(x))}
#'  @param nv the number of right singular vectors to be computed. Default is \code{min(dim(x))}
#'  @param k total number of components to return though the full variance (based on nu and nv) will still be returned (see \code{Dv.orig})
#'  @param tol A tolerance level for eliminating (tiny variance or negative or imaginary) components. Default is .Machine$double.eps
#'
#'  @return
#'  A list with seven elements:
#'  \item{u} Left singular vectors. A matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension c(n, nu) but also accounting for \code{tol}.
#'  \item{v} Right singular vectors. A matrix whose columns contain the left singular vectors of x, present if nv > 0. Dimension c(p, nv) but also accounting for \code{tol}.
#'  \item{p} Left generalized singular vectors. A vector containing the singular values of x of length min(n, p) but also accounting for \code{tol}.
#'  \item{q} Right generalized singular vectors. A vector containing the singular values of x of length min(n, p) but also accounting for \code{tol}.
#'  \item{d} a vector containing the singular values of x of length min(n, p) but also accounting for \code{tol}.
#'  \item{d.orig} a vector containing the singular values of x of length min(n, p) but also accounting for \code{tol}.
#'  \item{tau} a vector that contains the (original) explained variance per component.
#'
#'  @seealso \code{\link{tolerance.svd}} and \code{\link{svd}}
#'
#'  @examples
#'  ## an example with correspondence analysis.
#'  authors <- rbind(
#'    cbind(7836, 13112, 6026),
#'    cbind(53655, 102383, 42413),
#'    cbind(115615, 184541, 59226),
#'    cbind(161926, 340479, 62754),
#'    cbind(38177, 105101, 12670),
#'    cbind(46371, 58367, 14299)
#'    )
#'
#'  Observed <- authors/sum(authors)
#'  row.w <- rowSums(Observed)
#'    row.W <- diag(1/row.w)
#'  col.w <- colSums(Observed)
#'    col.W <- diag(1/col.w)
#'  Expected <- row.w %o% col.w
#'  Deviations <- Observed - Expected
#'  ca.res <- gsvd(Deviations,row.W,col.W)
#'    ## fi & fj in example deprecated because this is now part of the GSVD return.
#'  #fi <- row.W %*% ca.res$p %*% diag(ca.res$d)
#'  #fj <- col.W %*% ca.res$q %*% diag(ca.res$d)
#'
#'  @author Derek Beaton
#'  @keywords multivariate, diagonalization, eigen


gsvd <- function(DAT, LW=NaN, RW=NaN, nu= min(dim(DAT)), nv = min(dim(DAT)), k = 0, tol=.Machine$double.eps){

    ## probably need some rudimentary checks here.
    DAT <- as.matrix(DAT)

    ## I should check if LW and RW are rectangular. If they are, I should exit with an error.

    RW.is.vector <- LW.is.vector <- RW.is.nan <- LW.is.nan <- F
      ## clean this up to avoid the warnings...
    if( is.nan(LW) ){
      LW.is.nan <- T
    }else{
      if ( is.null(dim(LW)) & (length(LW) > 0) ) {
        LW.is.vector <- T
      }
      if(!LW.is.vector){
        if( isDiagonal.matrix(LW) ){
          LW <- diag(LW)
          LW.is.vector <- T
        }
      }
    }
    if( is.nan(RW) ){
      RW.is.nan <- T
    }else{
      if ( is.null(dim(RW)) & (length(RW) > 0) ) {
        RW.is.vector <- T
      }
      if(!RW.is.vector){
        if( isDiagonal.matrix(RW) ){
          RW <- diag(RW)
          RW.is.vector <- T
        }
      }
    }


    if( LW.is.vector ){
      DAT <- matrix(sqrt(LW),nrow=nrow(DAT),ncol=ncol(DAT),byrow=F) * DAT
    }else if(!LW.is.nan){
      DAT <- power.rebuild_matrix(LW, power = 1/2) %*% DAT
    }
    if( RW.is.vector ){
      DAT <- DAT * matrix(sqrt(RW),nrow=nrow(DAT),ncol=ncol(DAT),byrow=T)
    }else if(!RW.is.nan){
      DAT <- DAT %*% power.rebuild_matrix(RW, power = 1/2)
    }

  ## I also need to skip over this computation if LW or RW are either empty or all 1s
  #dat.for.svd <- power.rebuild_matrix(LW, power = 1/2) %*% DAT %*% power.rebuild_matrix(RW, power = 1/2)

  if(k<=0){
    k <- min(nrow(DAT),ncol(DAT))
  }
  res <- tolerance.svd(DAT,nu=nu,nv=nv,tol=tol)

  d.orig <- res$d
  tau <- d.orig^2/sum(d.orig^2)
  comp.ret <- min(length(d.orig),k)

  d <- d.orig[1:comp.ret]
          ## this should protect against the rank 1 where it's just a vector.
  res$u <- as.matrix(res$u[,1:comp.ret])
  res$v <- as.matrix(res$v[,1:comp.ret])


    ## I also need to skip over this computation if LW or RW are either empty or all 1s
  if(LW.is.vector){
    p <- matrix(1/sqrt(LW),nrow=nrow(res$u),ncol=ncol(res$u),byrow=F) * res$u
    fi <- matrix(LW,nrow=nrow(p),ncol=ncol(p),byrow=F) * p * matrix(d,nrow(p),ncol(p),byrow=T)
  }else if(!LW.is.nan){
    p <- power.rebuild_matrix(LW, power = -1/2) %*% res$u
    fi <- LW %*% p * matrix(d,nrow(p),ncol(p),byrow=T)
  }else{
    p <- res$u
    fi <- p * matrix(d,nrow(p),ncol(p),byrow=T)
  }

  if(RW.is.vector){
    q <- matrix(1/sqrt(RW),nrow=nrow(res$v),ncol=ncol(res$v),byrow=F) * res$v
    fj <- matrix(RW,nrow=nrow(q),ncol=ncol(q),byrow=F) * q * matrix(d,nrow(q),ncol(q),byrow=T)
  }else if(!RW.is.nan){
    q <- power.rebuild_matrix(RW, power = -1/2) %*% res$v
    fj <- RW %*% q * matrix(d,nrow(q),ncol(q),byrow=T)
  }else{
    q <- res$v
    fj <- q * matrix(d,nrow(q),ncol(q),byrow=T)
  }

  rownames(fi) <- rownames(res$u) <- rownames(p) <- rownames(DAT)
  rownames(fj) <- rownames(res$v) <- rownames(q) <- colnames(DAT)

  return(list(fi = fi, fj = fj, p = p, q = q, u = res$u, v = res$v, d = d, d.orig = d.orig, tau = tau))
}

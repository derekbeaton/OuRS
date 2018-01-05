#' @export
#'
#' @title \code{tolerance.svd}: automatically truncate spurious (tiny variance or negative or imaginary) components.
#'
#' @description \code{tolerance.svd} eliminates likely spurious components: any eigenvalue (squared singular value) below a tolerance level is elminated.
#'    The (likely) spurious singular values and vectors are then eliminated from \code{$u}, \code{$d}, and \code{$v}.
#'    Additionally, all values in \code{abs($u)} or \code{abs($v)} that fall below the \code{tol} are set to 0.
#'
#' @param x A data matrix of size for input to the singular value decomposition (\code{\link{svd}})
#' @param nu The number of left singular vectors to be computed. Default is \code{min(dim(x))}
#' @param nv The number of right singular vectors to be computed. Default is \code{min(dim(x))}
#' @param tol Default is \code{.Machine$double.eps}. A parameter with two roles: A tolerance level for 1: eliminating (tiny variance or negative or imaginary) components and 2: converting all values < tol to 0 in \code{u} and \code{v}
#'
#' @return A list with three elements (like \code{svd}):
#'  \item{d}{ A vector containing the singular values of x > \code{tol}.}
#'  \item{u}{ A matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension \code{min(c(nrow(x), nu, length(d))}.}
#'  \item{v}{ A matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension \code{min(c(ncol(x), nv, length(d))}.}
#'
#' @seealso \code{\link{svd}}
#'
#' @examples
#'  hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#'  X <- hilbert(9)[, 1:6]
#'  s_asis <- tolerance.svd(X)
#'  s_sqrt.Machine <- tolerance.svd(X,tol=sqrt(.Machine$double.eps))
#'  s_000001 <- tolerance.svd(X,tol=.000001)
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen

tolerance.svd <- function(x, nu=min(dim(x)), nv=min(dim(x)), tol=.Machine$double.eps) {	## consider increasing the tolerance.

    ##nu and nv are pass through values.
  svd.res <- svd(x, nu = nu, nv = nv)
  if(any(unlist(lapply(svd.res$d,is.complex)))){
    stop("tolerance.svd: Singular values ($d) are complex.")
  }
  svs.to.keep <- which(!(svd.res$d^2 < tol))
  svd.res$d <- svd.res$d[svs.to.keep]

  if(nu >= length(svs.to.keep)){
    svd.res$u <- as.matrix(svd.res$u[,svs.to.keep])
  }else{
    svd.res$u <- as.matrix(svd.res$u[,1:nu])
  }
  rownames(svd.res$u) <- rownames(x)
  svd.res$u[ abs(svd.res$u) < tol ] <- 0


  if(nv >= length(svs.to.keep)){
    svd.res$v <- as.matrix(svd.res$v[,svs.to.keep])
  }else{
    svd.res$v <- as.matrix(svd.res$v[,1:nv])
  }
  rownames(svd.res$v) <- colnames(x)
  svd.res$v[ abs(svd.res$v) < tol ] <- 0


  return(svd.res)
}

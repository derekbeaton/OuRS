### NOTE TO SELF: Should I eliminate complex eigenvalues (and vectors)?
### FOR NOW: STOP if complex SVs/Eigens are encountered.

### ANOTHER NOTE: 
	## Maybe this function should have k= and return the low rank parts instead of doing that elsewhere.

#  "Tolerance" SVD: A SVD function that automatically eliminates (likely) spurious components/sources of variance
#    These eliminated components are those that fall below a specified tolerance, currently defaulted to .Machine$double.eps
#'
#'  @export
#'
#'  @title \code{tolerance.svd}: automatically truncate spurious (tiny variance or negative or imaginary) components.
#'
#'  @description \code{tolerance.svd} eliminates likely spurious components: any eigenvalue (squared singular value) below a tolerance level is elminated.
#'    The (likely) spurious singular values and vectors are then eliminated from \code{$u}, \code{$d}, and \code{$v}.
#'    Additionally, all values in \code{abs($u)} or \code{abs($v)} that fall below the \code{tol} are set to 0.
#'
#'  @param x A data matrix for input to the singular value decomposition (\code{\link{svd}})
#'  @param nu the number of left singular vectors to be computed. Default is \code{min(dim(x))}
#'  @param nv the number of right singular vectors to be computed. Default is \code{min(dim(x))}
#'  @param tol A tolerance level for eliminating (tiny variance or negative or imaginary) components. Default is .Machine$double.eps
#'
#'  @return
#'  A list with three elements (like \code{svd}):
#'  \item{d} a vector containing the singular values of x of length min(n, p) but also accounting for \code{tol}.
#'  \item{u} a matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension c(n, nu) but also accounting for \code{tol}.
#'  \item{v} a matrix whose columns contain the left singular vectors of x, present if nv > 0. Dimension c(p, nv) but also accounting for \code{tol}.
#'
#'  @seealso \code{\link{svd}}
#'
#'  @examples
#'  hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#'  X <- hilbert(9)[, 1:6]
#'  s_asis <- tolerance.svd(X)
#'  s_sqrt.Machine <- tolerance.svd(X,tol=sqrt(.Machine$double.eps))
#'  s_000001 <- tolerance.svd(X,tol=.000001)
#'
#'  @author Derek Beaton
#'  @keywords multivariate, diagonalization, eigen

tolerance.svd <- function(x, nu=min(dim(x)), nv=min(dim(x)), tol=.Machine$double.eps) {	## consider increasing the tolerance.

    ##nu and nv are pass through values.
  svd.res <- svd(x, nu = nu, nv = nv)
  if(any(unlist(lapply(svd.res$d,is.complex)))){
    stop("tolerance.svd: Singular values ($d) are complex.")
  }
  comps.to.keep <- which(!(svd.res$d^2 < tol))

  svd.res$u <- as.matrix(svd.res$u[,comps.to.keep])
  rownames(svd.res$u) <- rownames(x)
  svd.res$v <- as.matrix(svd.res$v[,comps.to.keep])
  rownames(svd.res$v) <- colnames(x)
  svd.res$d <- svd.res$d[comps.to.keep]

  svd.res$u[ abs(svd.res$u) < tol ] <- 0
  svd.res$v[ abs(svd.res$v) < tol ] <- 0

  return(svd.res)

}

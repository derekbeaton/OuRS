#' @export
#'
#' @title Principal components analysis
#'
#' @description
#' \code{pca} performs principal components analysis of a data matrix \code{DATA}.
#'
#' @param DATA a data matrix to decompose
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param k total number of components to return (see \code{\link[GSVD]{gsvd}}).
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#' @return A list with eleven elements by way of \code{\link[GSVD]{gsvd}}:
#' \item{d.orig}{A vector containing the singular values of X above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of X above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}).}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values of X}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values of X}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DATA)} by k.}
#' \item{p}{Left (rows) generalized singular vectors. Dimensions are \code{nrow(DATA)} by k.}
#' \item{fi}{Left (rows) component scores. Dimensions are \code{nrow(DATA)} by k.}
#' \item{v}{Right (columns) singular vectors. Dimensions are \code{ncol(DATA)} by k.}
#' \item{q}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(DATA)} by k.}
#' \item{fj}{Right (columns) component scores. Dimensions are \code{ncol(DATA)} by k.}
#'
#' @seealso \code{\link{gsvd}}
#'
#'
#' @author Derek Beaton


pca <- function(DATA, center = T, scale = T, k = 0, tol = sqrt(.Machine$double.eps)){

  res <- gsvd(ours_scale(DATA, center = center, scale = scale), k = k, tol = tol)
  res

}



#' @export
#'
#' @title Correspondence analysis
#'
#' @description
#' \code{ca} performs correspondence analysis of a data matrix \code{DATA}.
#'
#' @param DATA a data matrix to decompose
#' @param k total number of components to return (see \code{\link{gsvd}}).
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#' @return A list with eleven elements by way of \code{\link[GSVD]{gsvd}}:
#' \item{d.orig}{A vector containing the singular values of X above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of X above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}).}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values of X}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values of X}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DATA)} by k.}
#' \item{p}{Left (rows) generalized singular vectors. Dimensions are \code{nrow(DATA)} by k.}
#' \item{fi}{Left (rows) component scores. Dimensions are \code{nrow(DATA)} by k.}
#' \item{v}{Right (columns) singular vectors. Dimensions are \code{ncol(DATA)} by k.}
#' \item{q}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(DATA)} by k.}
#' \item{fj}{Right (columns) component scores. Dimensions are \code{ncol(DATA)} by k.}
#'
#' @seealso \code{\link{gsvd}}
#'
#'
#' @author Derek Beaton


ca <- function(DATA, k = 0, tol = sqrt(.Machine$double.eps)){

  sum.data <- sum(DATA)
  wi <- rowSums(DATA)/sum.data
  wj <- colSums(DATA)/sum.data

  res <- gsvd( ((DATA/sum.data) - (wi %o% wj)), 1/wi, 1/wj, k = k, tol = tol )

  res
}

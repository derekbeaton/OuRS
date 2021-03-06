% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generalized_mcd.R
\name{generalized_mcd}
\alias{generalized_mcd}
\title{Generalized minimum covariance determinant (GMCD) for non-continuous data}
\usage{
generalized_mcd(
  DATA,
  alpha = 0.75,
  num.subsets = 500,
  max.total.iters = num.subsets * 20,
  top.sets.percent = 0.05,
  tol = .Machine$double.eps
)
}
\arguments{
\item{DATA}{a data matrix (of presumably all transformed data)}

\item{alpha}{numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point}

\item{num.subsets}{numeric. The number of initial subsamples to start the MCD algorithm with}

\item{max.total.iters}{numeric. The total number of iterations allowed for the MCD search}

\item{top.sets.percent}{numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)}

\item{tol}{default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).}
}
\value{
The 'OuRS MCD' object: a list of three lists:
\strong{cov} a list for the robust covariance structure items
\itemize{
  \item loadings - a matrix of loadings via CA from the robust covariance matrix
  \item singular.values - a vector of singular values via CA from the robust covariance matrix
}
\strong{dists}
\itemize{
  \item robust_mahal_dists - Robust Mahalanobis distances from the robust covariance matrix
  \item robust_score_dists - Robust score distances (computed from component scores) from the robust covariance matrix
  \item mahal_dists - Mahalanobis distances
  \item{score_dists} {Score distances (computed from component scores)}
}
\strong{det.samps}
\itemize{
  \item dets - A numeric vector. The \code{top.sets.percent} determinants in ascending order (from minimum determinant upwards) that reflects the \code{top.sets.percent} best determinants from the MCD search
  \item samples - A numeric matrix. The \code{top.sets.percent} subsamples in to compute the determinants (in \code{dets})
}
}
\description{
\code{generalized_mcd} performs the GMCD of a data matrix \code{DATA}.
}
\details{
the input parameter \code{DATA} are assumed to be transformed into disjunctive data (see \code{\link{disjunctive_coding}}) or an analog for continuous and/or ordinal data (see \code{\link{escofier_coding}}, \code{\link{thermometer_coding}}, \code{\link{mixed_data_coding}}).
Generally, \code{DATA} should have the same properties as a data matrix that would be analyzed by Correspondence Analysis.
}
\seealso{
\code{\link{categorical_mcd}}, \code{\link{ordinal_mcd}}, \code{\link{mixed_data_mcd}}, \code{\link{continuous_mcd}}
}
\author{
Derek Beaton
}

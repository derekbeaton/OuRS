#' @title Minimum covariance determinant (MCD) for continuous data
#'
#' @description
#' \code{continuous_mcd} performs the MCD of a data matrix \code{DATA}.
#' 
#'
#' @param DATA a data matrix (of presumably all continuous data)
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param allow_collinearity logical (boolean). Default is \code{FALSE} which does not allow a matrix to be singular (i.e., all eigenvalues must be > 0).
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the MCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the MCD search
#' @param top.sets.percent numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#' @return The 'OuRS MCD' object: a list of three lists:
#' \item \strong{cov:} a list for the robust covariance structure items
#' \itemize{
#'   \item{loadings:} {a matrix of loadings via PCA from the robust covariance matrix}
#'   \item{singular.values:} {a vector of singular values via PCA from the robust covariance matrix}
#'   \item{center:} {a vector of the column-wise centers for the final subsample that produces the robust covariance matrix}
#'   \item{scale:} {a vector of the column-wise scale for the final subsample that produces the robust covariance matrix}
#' }
#' \item \strong{dists}
#' \itemize{
#'   \item{robust_mahal_dists} {Robust Mahalanobis distances from the robust covariance matrix}
#'   \item{robust_score_dists} {Robust score distances (computed from component scores) from the robust covariance matrix}
#'   \item{mahal_dists} {Mahalanobis distances}
#'   \item{score_dists} {Score distances (computed from component scores)}
#' }
#' \item \strong{det.samps}
#' \itemize{
#'   \item{dets:} {A numeric vector. The \code{top.sets.percent} determinants in ascending order (from minimum determinant upwards) that reflects the \code{top.sets.percent} best determinants from the MCD search}
#'   \item{samples:} {A numeric matrix. The \code{top.sets.percent} subsamples in to compute the determinants (in \code{dets})}
#' }
#'
#' @seealso \code{\link{categorical_mcd}}, \code{\link{ordinal_mcd}}, \code{\link{mixed_data_mcd}}, \code{\link{generalized_mcd}}
#'
#' @examples
#'
#' @author Derek Beaton
#' @export

continuous_mcd <- function(DATA, center=T, scale=F, allow_collinearity=F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  ## data size check
  if(ncol(DATA) > (nrow(DATA)*.9)){
  	stop("continuous_mcd: the column-to-row ratio is too high so 'mcd' cannot be performed")
  }
  
  ## stop if anything is NA; ask that they handle NAs outside of here
  if( any(is.na(DATA)) | any(is.infinite(DATA)) | any(is.nan(DATA)) | any(is.null(DATA)) ){
    stop("continuous_mcd: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }
  
  ## collinearity check
    ## NOTE: the DetMCD.m does exactly this via classSVD.m
  svd.res <- tolerance_svd(ours_scale(DATA, center=center, scale=scale), tol = NA)
  if( any(svd.res$d^2 < tol) ) {
    if(!allow_collinearity){
      stop("continuous_mcd: Matrix is singular. It is likely that data are collinear where some variables are combinations of other variables.")
    }
  }
  
  ## sample finder
  mcd.samples <- continuous_mcd_search_for_sample(DATA, center=center, scale=scale, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent,tol=tol)

  ## get sample distances
  tsvd.res <- tolerance_svd(ours_scale(DATA, center=center, scale=scale), tol = tol)
  mahals <- rowSums(tsvd.res$u^2)
  chis <- rowSums( (tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T))^2 )

  ## get robust mean & cov (loadings)
  rob.sample <- ours_scale(DATA[mcd.samples$final_subsamples[1,],],center=center,scale=scale)
  rob.center <- attributes(rob.sample)$`scaled:center`
  rob.scale <- attributes(rob.sample)$`scaled:scale`
  robust.tsvd.res <- tolerance_svd(rob.sample,tol = tol)

  ## get robust distances
  robust.vectors.and.scores <- continuous_scores_dists(DATA, center=rob.center, scale=rob.scale,robust.tsvd.res$v,robust.tsvd.res$d)


  ## I should also return the actual scores, too
  res <- list(
    cov = list(loadings = robust.tsvd.res$v,
               singular.values = robust.tsvd.res$d,
               center = rob.center,
               scale = rob.scale
    ),
    dists = list(robust_mahal_dists = robust.vectors.and.scores$projected_mahal_dists,
                 rob_score_dists = robust.vectors.and.scores$projected_score_dists,
                 mahal_dists = mahals,
                 score_dists = chis),
    det.samps = list(dets = mcd.samples$final_determinants,
                     samples = mcd.samples$final_subsamples)
  )
  class(res) <- c("list", "OuRS", "MCD", "continuous")

  return(res)
}


### eventually, these things should be class based
#### I have accepted the fact that all of this will be completely re-written some day
##### perhaps for ExPo2?

#' @title Compute scores and distances with projection
#'
#' @description Computes projected: singular vectors, Mahalanobis distance, component scores, and score distances. 
#' 
#' @details This approach uses loadings (singular vectors) and singular values to compute scores and distances for \code{DATA}. For use with the MCD, \code{loadings} and \code{singular.values} come from a robust covariance matrix.
#'
#' @param DATA a data matrix (of presumably all continuous data)
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param loadings a
#' @param singular.values a
#'
#' @return a list with four items. All items are for the rows of \code{DATA} and computed through projection (via \code{loadings} and \code{singular.vectors})
#' \item{projected_u:} {Projected singular vectors}
#' \item{projected_fi:} {Projected component scores}
#' \item{projected_mahal_dists:} {Mahalanobis distances (computed as \code{rowSums(projected_u^2)})}
#' \item{projected_score_dists:} {Score distances (computed as \code{rowSums(projected_fi^2)})}
#'
#' @seealso \code{\link{continuous_mcd}} and \code{\link{generalized_scores_dists}}
#'
#' @examples
#'
#' @author Derek Beaton
#' @export

continuous_scores_dists <- function(DATA, center=T, scale=F, loadings, singular.values){

  projected_fi <- (ours_scale(DATA,center=center,scale=scale) %*% loadings)
  projected_u <- sweep(projected_fi,2,singular.values,"/")

  projected_mahal_dists <- rowSums(projected_u^2)
  projected_score_dists <- rowSums(projected_fi^2)

  return( list(projected_mahal_dists = projected_mahal_dists, projected_score_dists = projected_score_dists, projected_u = projected_u, projected_fi = projected_fi) )

}


#' @title Minimum covariance determinant (MCD) search for subsample in continuous data
#'
#' @description Performs the search for the subsample for MCD of a data matrix \code{DATA}.
#'
#' @param DATA a data matrix (of presumably all continuous data)
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the MCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the MCD search
#' @param top.sets.percent numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#' @return a list with a vector and a matrix
#' \item{final_determinants: }{a vector of length \code{round(num.subsets * top.sets.percent)} with determinants produced from the search. Listed in ascending order, so the first element is the minimum determinant}
#' \item{final_subsamples: }{a matrix \code{round(num.subsets * top.sets.percent)} rows with subsamples produced from the search that produce the \code{final_determinants}. Listed in descending order, so the first row is the subsample that produces the minimum determinant}
#'
#'
#' @author Derek Beaton
#' 
#' @noRd


continuous_mcd_search_for_sample <- function(DATA, center=T, scale=F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  # h.alpha.n now handles changing of alpha values outside of the range of [.5,1]
  h.size <- h.alpha.n(alpha, nrow(DATA), ncol(DATA))
  max.det.iters <- round(max.total.iters / num.subsets)
  
  dets <- vector("numeric", num.subsets)
  orders <- matrix(NA,num.subsets,h.size)
  for(i in 1:num.subsets){
    
    findInit <- T
    init.size <- min(dim(DATA))+1
    
    ## this step isn't necessary, it's just in the spirit of the MCD algorithm
    ## we start with the minimum possible set size and work from there
    ## we could just dive into the c_step here (but that's for later)
    while( findInit ){
      
      init.samp <- sort(sample(nrow(DATA),init.size))
      init.norm <- ours_scale(DATA[init.samp,],center,scale)
      init.svd <- tolerance_svd(init.norm, tol = tol)
      init.mds <- round(rowSums(init.svd$u^2),digits=8)
      
      if(length(unique(init.mds)) < 2){
        init.size <- init.size + 1
      }else{
        sup.scores.and.dists <- continuous_scores_dists(DATA, attributes(init.norm)$`scaled:center`, attributes(init.norm)$`scaled:scale`, init.svd$v, init.svd$d)
        samp.config <- sort(order(sup.scores.and.dists$projected_mahal_dists)[1:h.size])
        findInit <- F
      }
    }
    
    min.info <- continuous_c_step(DATA, center=center, scale=scale, observations_subsample = samp.config, max.iters = max.det.iters)
    dets[i] <- min.info$minimum_determinant
    orders[i,] <- min.info$observations_subsample
  }
  
  perc.cut <- round(num.subsets * top.sets.percent)
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique(unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),])
  
  final_determinants <- vector("numeric", nrow(final.configs))
  final_subsamples <- matrix(NA,nrow(final.configs),h.size)
  for( i in 1:nrow(final.configs)){
    
    min.info <- continuous_c_step(DATA, center = center, scale = scale, observations_subsample = final.configs[i,])
    final_determinants[i] <- min.info$minimum_determinant
    final_subsamples[i,] <- min.info$observations_subsample
    
  }
  
  best.order <- order(final_determinants)
  final_determinants <- final_determinants[best.order]
  final_subsamples <- final_subsamples[best.order,]
  
  return( list(final_determinants = final_determinants, final_subsamples = final_subsamples) )
  
}


#' @title C step of the MCD
#' 
#' @description Performs the C ("concentration") step in the search for a sample of observations that produces a minimum determinant from a covariance matrix
#' 
#' @param DATA a data matrix (of presumably all continuous data)
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param observations_subsample a numeric vector. Indicates which observations to use as the initial subsample in the C-step
#' @param max.iters numeric. Default is 100. Indicates how many iterations of the C-step to perform before stopping
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#' @return a list with two items
#' \item{observations_subsample: }{a numeric vector to for the subsample that produces a minimum determinant}
#' \item{minimum_determinant: }{the determinant of the covariance matrix produced by \code{observations_subsample}}
#'
#' @author Derek Beaton
#'
#' @noRd

continuous_c_step <- function(DATA, center=T, scale=F, observations_subsample, max.iters=100, tol=.Machine$double.eps){

  old_determinant <- Inf
  old_center <- NaN
  old_loadings <- matrix(NaN,0,0)
  new_subsample <- old_subsample <- observations_subsample

  for(i in 1:max.iters){

    sub.data <- DATA[new_subsample,]
    sub.data.normed <- ours_scale(sub.data,center=center,scale=scale)
    svd.res <- tolerance_svd(sub.data.normed, tol = tol)

    new_center <- attributes(sub.data.normed)$`scaled:center`
    new_scale <- attributes(sub.data.normed)$`scaled:scale`
    # new_determinant <- geometric_mean(svd.res$d^2)
    new_determinant <- exp(mean(log(svd.res$d^2)))

    if( (new_determinant <= old_determinant) & (!isTRUE(all.equal(new_determinant,0,tolerance=tol))) ){
      if( center.sigma_checker(old_center, new_center, old_loadings, svd.res$v, tol=tol) & isTRUE(all.equal(new_determinant, old_determinant, tolerance= tol)) ){
        return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
      }
      else{

        old_determinant <- new_determinant
        old_center <- new_center
        old_loadings <- svd.res$v
        old_subsample <- new_subsample

        sup.scores.and.dists <- continuous_scores_dists(DATA, new_center, new_scale, svd.res$v, svd.res$d)
        new_subsample <- sort(order(sup.scores.and.dists$projected_mahal_dists)[1:length(observations_subsample)])

        if( isTRUE(all.equal(sort(new_subsample),sort(old_subsample))) ){
          return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
        }
      }
    }else{
      return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
    }
  }

  return( list(observations_subsample = new_subsample, minimum_determinant = new_determinant) )
}

### a potential speed up and big change:
  ### from the current usage, I'm not allowing external/reference weights or expected matrices
  ### everything here is computed directly from the data
  ### so I should be able to overhaul all the cals that take in the wZx, m, w, and such, and thenn only pass in DATA
    #### I do suppose that it's possibly faster to pass that information in, but then we have duplicate information all over the place
    #### it does seem like I only need the ca_preproc info one time...
    #### so I should make a choice: simplify here and force the use of these in one way or
    #### allow for the more generalized approach, which requires more parameters to pass, but only a singular computation of ca_preproc

### quick exec. decision: it all stays for now because it is inevitable that there will be a rewrite for optimization


### the categorical version should explicitly handle the disjunctive transform
## same for ordinal and mixed data
## but all of these functions just pass the transformed data off to a generalized_mcd() which just takes data
  ## and it performs it as defined in the paper.


#' @title Generalized minimum covariance determinant (GMCD) for ordinal data
#'
#' @description
#' \code{ordinal_mcd} performs the GMCD of a data matrix \code{DATA}.
#' 
#' @param DATA a numeric data matrix (of presumably all ordinal data)
#' @param mins a vector to denote the expected minimum values for \code{DATA}. If not defined, the observed minimums will be used
#' @param maxs a vector to denote the expected maximum values for \code{DATA}. If not defined, the observed maximums will be used
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the MCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the MCD search
#' @param top.sets.percent numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#'
#' @return The 'OuRS MCD' object: a list of three lists:
#' \strong{cov} a list for the robust covariance structure items
#' \itemize{
#'   \item loadings - a matrix of loadings via CA from the robust covariance matrix
#'   \item singular.values - a vector of singular values via CA from the robust covariance matrix
#' }
#' \strong{dists}
#' \itemize{
#'   \item robust_mahal_dists - Robust Mahalanobis distances from the robust covariance matrix
#'   \item robust_score_dists - Robust score distances (computed from component scores) from the robust covariance matrix
#'   \item mahal_dists - Mahalanobis distances
#'   \item score_dists - Score distances (computed from component scores)
#' }
#' \strong{det.samps}
#' \itemize{
#'   \item dets - A numeric vector. The \code{top.sets.percent} determinants in ascending order (from minimum determinant upwards) that reflects the \code{top.sets.percent} best determinants from the MCD search
#'   \item samples - A numeric matrix. The \code{top.sets.percent} subsamples in to compute the determinants (in \code{dets})
#' }
#'
#' @seealso \code{\link{continuous_mcd}}, \code{\link{mixed_data_mcd}}, \code{\link{categorical_mcd}}, \code{\link{generalized_mcd}}
#'
#'
#' @author Derek Beaton
#' @export
#' 
ordinal_mcd <- function(DATA, mins = NULL, maxs = NULL, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  # if(make_data_doubled){
    DATA <- thermometer_coding(DATA, mins = mins, maxs = maxs)
  # }

  
  res <- generalized_mcd(DATA, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent, tol=tol)
  class(res) <- c("list", "OuRS", "MCD", "generalized", "ordinal")
  return(res)

}


## data types only allow for n, c, z, o, and x
#' @title Generalized minimum covariance determinant (GMCD) for mixed data
#'
#' @description
#' \code{mixed_data_mcd} performs the GMCD of a data matrix \code{DATA}.
#' 
#' @details the input parameter \code{DATA} must be a \code{data.frame}, where ordinal and continuous variables are \code{numeric} and categorical variables are \code{character}. Any variables denoted with "x" (i.e., "do nothing") must also be \code{numeric}.
#' See \code{\link{mixed_data_coding}}
#'
#' @param DATA a data frame with character columns (categorical) and numeric otherwise
#' @param column.types a vector to denote the data type of each column of \code{DATA}. Options are "n" for categorical (nominal), "c" for continuous (centering only), "z" for continuous (centering and scaling), "o" for ordinal, and "x" which does nothing to those respective columns
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the MCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the MCD search
#' @param top.sets.percent numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#'
#' @return The 'OuRS MCD' object: a list of three lists:
#' \strong{cov} a list for the robust covariance structure items
#' \itemize{
#'   \item loadings - a matrix of loadings via CA from the robust covariance matrix
#'   \item singular.values - a vector of singular values via CA from the robust covariance matrix
#' }
#' \strong{dists}
#' \itemize{
#'   \item robust_mahal_dists - Robust Mahalanobis distances from the robust covariance matrix
#'   \item robust_score_dists - Robust score distances (computed from component scores) from the robust covariance matrix
#'   \item mahal_dists - Mahalanobis distances
#'   \item score_dists - Score distances (computed from component scores)
#' }
#' \strong{det.samps}
#' \itemize{
#'   \item dets - A numeric vector. The \code{top.sets.percent} determinants in ascending order (from minimum determinant upwards) that reflects the \code{top.sets.percent} best determinants from the MCD search
#'   \item samples - A numeric matrix. The \code{top.sets.percent} subsamples in to compute the determinants (in \code{dets})
#' }
#'
#' @seealso \code{\link{continuous_mcd}}, \code{\link{ordinal_mcd}}, \code{\link{categorical_mcd}}, \code{\link{generalized_mcd}}
#'
#'
#' @author Derek Beaton
#' @export
#' 
mixed_data_mcd <- function(DATA, column.types=rep("x",ncol(DATA)), alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  
  DATA <- mixed_data_coding(DATA, column.types)
  res <- generalized_mcd(DATA, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent, tol=tol)
  class(res) <- c("list", "OuRS", "MCD", "generalized", "mixed")
  return(res)
  
}


#' @title Generalized minimum covariance determinant (GMCD) for categorical data
#'
#' @description
#' \code{categorical_mcd} performs the GMCD of a data matrix \code{DATA}.
#' 
#' @details the input parameter \code{DATA} must be a matrix of categorical data. Each level for each column (variable) will be represented through disjunctive coding. See \code{\link{disjunctive_coding}}
#'
#' @param DATA a data matrix (of presumably all categorical data)
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the MCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the MCD search
#' @param top.sets.percent numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#'
#' @return The 'OuRS MCD' object: a list of three lists:
#' \strong{cov} a list for the robust covariance structure items
#' \itemize{
#'   \item loadings - a matrix of loadings via CA from the robust covariance matrix
#'   \item singular.values - a vector of singular values via CA from the robust covariance matrix
#' }
#' \strong{dists}
#' \itemize{
#'   \item robust_mahal_dists - Robust Mahalanobis distances from the robust covariance matrix
#'   \item robust_score_dists - Robust score distances (computed from component scores) from the robust covariance matrix
#'   \item mahal_dists - Mahalanobis distances
#'   \item score_dists - Score distances (computed from component scores)
#' }
#' \strong{det.samps}
#' \itemize{
#'   \item dets - A numeric vector. The \code{top.sets.percent} determinants in ascending order (from minimum determinant upwards) that reflects the \code{top.sets.percent} best determinants from the MCD search
#'   \item samples - A numeric matrix. The \code{top.sets.percent} subsamples in to compute the determinants (in \code{dets})
#' }
#'
#' @seealso \code{\link{continuous_mcd}}, \code{\link{ordinal_mcd}}, \code{\link{mixed_data_mcd}}, \code{\link{generalized_mcd}}
#'
#'
#' @author Derek Beaton
#' @export
#' 
categorical_mcd <- function(DATA, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  # if(make_data_disjunctive){
    DATA <- disjunctive_coding(DATA)
  # }

  res <- generalized_mcd(DATA, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent, tol=tol)
  class(res) <- c("list", "OuRS", "MCD", "generalized", "categorical")
  return(res)
  
}


### techinically, the generalized MCD is CA-based, so it has to assume we have data that are CA-friendly
#' @title Generalized minimum covariance determinant (GMCD) for non-continuous data
#'
#' @description
#' \code{generalized_mcd} performs the GMCD of a data matrix \code{DATA}.
#' 
#' @details the input parameter \code{DATA} are assumed to be transformed into disjunctive data (see \code{\link{disjunctive_coding}}) or an analog for continuous and/or ordinal data (see \code{\link{escofier_coding}}, \code{\link{thermometer_coding}}, \code{\link{mixed_data_coding}}).
#' Generally, \code{DATA} should have the same properties as a data matrix that would be analyzed by Correspondence Analysis.
#'
#' @param DATA a data matrix (of presumably all transformed data)
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the MCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the MCD search
#' @param top.sets.percent numeric. A value within (0,1] for the number of samples and determinants to return. Returned results are in ascending order of determinants (minimum is the first element)
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components (see \code{\link{gsvd}}).
#'
#'
#' @return The 'OuRS MCD' object: a list of three lists:
#' \strong{cov} a list for the robust covariance structure items
#' \itemize{
#'   \item loadings - a matrix of loadings via CA from the robust covariance matrix
#'   \item singular.values - a vector of singular values via CA from the robust covariance matrix
#' }
#' \strong{dists}
#' \itemize{
#'   \item robust_mahal_dists - Robust Mahalanobis distances from the robust covariance matrix
#'   \item robust_score_dists - Robust score distances (computed from component scores) from the robust covariance matrix
#'   \item mahal_dists - Mahalanobis distances
#'   \item{score_dists} {Score distances (computed from component scores)}
#' }
#' \strong{det.samps}
#' \itemize{
#'   \item dets - A numeric vector. The \code{top.sets.percent} determinants in ascending order (from minimum determinant upwards) that reflects the \code{top.sets.percent} best determinants from the MCD search
#'   \item samples - A numeric matrix. The \code{top.sets.percent} subsamples in to compute the determinants (in \code{dets})
#' }
#'
#' @seealso \code{\link{categorical_mcd}}, \code{\link{ordinal_mcd}}, \code{\link{mixed_data_mcd}}, \code{\link{continuous_mcd}}
#'
#'
#' @author Derek Beaton
#' @export
#' 

### I have to change all the chis to score_distance and md to mahal or something
generalized_mcd <- function(DATA, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  if(ncol(DATA) > (nrow(DATA)*.9)){
    stop("generalized_mcd: the column-to-row ratio is too high so 'mcd' cannot be performed")
  }
  
  ## stop if anything is NA; ask that they handle NAs outside of here
  if( any(is.na(DATA)) | any(is.infinite(DATA)) | any(is.nan(DATA)) | any(is.null(DATA)) ){
    stop("generalized_mcd: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }
  
  ## sample finder
  mcd.samples <- generalized_mcd_search_for_sample(DATA, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent)
  best.sample <- mcd.samples$final_subsamples[1,]

  preproc.DATA <- ca_preproc(DATA, compact = T) ## this could be more efficient...

  ca.res <- ca(DATA)
  mahals <- rowSums(ca.res$u^2)
  chis <- rowSums(ca.res$fi^2)


  ## get robust mean & cov (loadings)
  rob.sample <- preproc.DATA$weightedZx[best.sample,]
  robust.tsvd.res <- tolerance_svd(rob.sample,tol=tol)

  ## call to function that computes robust mahal
  ## I may not need to pass all of this in
  robust.vectors.and.scores <- generalized_scores_dists(DATA, preproc.DATA$m, preproc.DATA$w, robust.tsvd.res$v, robust.tsvd.res$d)

  res <- list(
    cov = list(loadings = robust.tsvd.res$v,
               singular.values = robust.tsvd.res$d
               # row_weights = preproc.DATA$m,
               # column_weights = preproc.DATA$w
    ),
    dists = list(robust_mahal_dists = robust.vectors.and.scores$projected_mahal_dists,
                 rob_score_dists = robust.vectors.and.scores$projected_score_dists,
                 mahal_dists = mahals,
                 score_dists = chis),
    det.samps = list(dets = mcd.samples$final_determinants,
                     samples = mcd.samples$final_subsamples)
  )


  class(res) <- c("list", "OuRS", "MCD", "generalized")
  return(res)

}

## same as the above two: this needs to be the generalized version, it doesn't care what data they were
  # formerly cat.sup.fi.u
#' @title Compute scores and distances with projection for generalized MCD
#'
#' @description Computes projected: singular vectors, Mahalanobis distance, component scores, and score distances. 
#' 
#' @details This approach uses loadings (singular vectors) and singular values to compute scores and distances for \code{DATA}. For use with the GMCD, \code{loadings} and \code{singular.values} come from a robust covariance matrix, where \code{row.weights} and \code{col.weights} are the weights from preprocessing for Correspondence Analysis.
#' 
#'
#' @param DATA a data matrix (of presumably all continuous data)
#' @param row.weights a numeric vector of row weights (as in correspondence analysis)
#' @param col.weights a numeric vector of column weights (as in correspondence analysis)
#' @param loadings a numeric matrix that contains the loadings (singular vectors or eigenvectors of variables)  from a decomposed covariance matrix
#' @param singular.values a numeric vector that contains the singular values from a decomposed covariance matrix
#'
#' @return a list with four items. All items are for the rows of \code{DATA} and computed through projection (via \code{loadings} and \code{singular.vectors})
#' \itemize{
#' \item projected_u - Projected singular vectors
#' \item projected_fi - Projected component scores
#' \item projected_mahal_dists - Mahalanobis distances (computed as \code{rowSums(projected_u^2)})
#' \item projected_score_dists - Score distances (computed as \code{rowSums(projected_fi^2)})
#' }
#' @seealso \code{\link{generalized_mcd}} and \code{\link{continuous_scores_dists}}
#'
#' @author Derek Beaton
#' @export

generalized_scores_dists <- function(DATA, row.weights, col.weights, loadings, singular.values){

  ## NOT EFFICIENT. MAKE MORE EFFICIENT
  profiles <- DATA / rowSums(DATA)
  projected_fi <- profiles %*% sweep(loadings,1,sqrt(col.weights)/col.weights,"*")
  projected_u <- sweep(sweep(projected_fi,2,singular.values,"/"),1,sqrt(row.weights),"*")

  projected_mahal_dists <- rowSums(projected_u^2)
  projected_score_dists <- rowSums(projected_fi^2)

  return( list(projected_mahal_dists = projected_mahal_dists, projected_score_dists = projected_score_dists, projected_u = projected_u, projected_fi = projected_fi) )

}




### techinically, the generalized MCD is CA-based, so it has to assume we have data that are CA-friendly
#### so this has the same requirement as generalized_mcd and works under the assumption that the data have been transformed appropriately

# formerly cat.mcd.find.sample
#' @title Generalized minimum covariance determinant (GMCD) search for subsample in non-continuous data
#'
#' @description Performs the search for the subsample for GMCD of a data matrix \code{DATA}.
#'
#' @param DATA a data matrix (of numeric but presumably non-continuous data that have been transformed)
#' @param alpha numeric. A value between .5 and 1 to select the size of the subsample based on a breakdown point
#' @param num.subsets numeric. The number of initial subsamples to start the GMCD algorithm with
#' @param max.total.iters numeric. The total number of iterations allowed for the GMCD search
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

generalized_mcd_search_for_sample <- function(DATA, alpha=.75,num.subsets=500,max.total.iters=num.subsets*20,top.sets.percent=.05,tol=.Machine$double.eps){
  
  ## I may not need this here.
    ### actually it is needed here.
    ### but is conceptually a place holder for the ability to use more flexible weights/reference structures
  preproc.DATA <- ca_preproc(DATA, compact = T)
  
  h.size <- h.alpha.n(alpha,nrow(DATA),ncol(DATA))
  max.det.iters <- round(max.total.iters / num.subsets)
  
  dets <- vector("numeric", num.subsets)
  orders <- matrix(NA,num.subsets,h.size)
  for(i in 1:num.subsets){
    
    findInit <- T
    init.size <- min(dim(DATA))+1
    
    while( findInit ){
      
      init.samp <- sort(sample(nrow(DATA),init.size))
      init.svd <- tolerance_svd(preproc.DATA$weightedZx[init.samp,],tol=tol)
      init.mds <- round(rowSums(init.svd$u^2),digits=8)	## do I need to round? ### I should probably use a tol parameter here...
      
      if(length(unique(init.mds)) < 2){
        init.size <- init.size + 1
      }else{
        
        ### generalized_scores_dists may not need all of these to be passed in...
        sup.scores.and.dists <- generalized_scores_dists(DATA, preproc.DATA$m, preproc.DATA$w, init.svd$v, init.svd$d)
        samp.config <- sort(order(sup.scores.and.dists$projected_mahal_dists)[1:h.size])
        findInit <- F
      }
    }
    
    ### generalized_c_step may not need all of these to be passed in...
    min.info <- generalized_c_step(DATA, weighted.deviations = preproc.DATA$weightedZx, row.weights = preproc.DATA$m, col.weights = preproc.DATA$w, observations_subsample = samp.config, max.iters = max.det.iters)
    dets[i] <- min.info$minimum_determinant
    orders[i,] <- min.info$observations_subsample
  }
  
  perc.cut <- round(num.subsets * top.sets.percent)
  unique.min.configs <- unique(orders[order(dets),])
  final.configs <- unique(unique.min.configs[1:min(nrow(unique.min.configs), perc.cut),])
  
  final_determinants <- vector("numeric", nrow(final.configs))
  final_subsamples <- matrix(NA,nrow(final.configs),h.size)
  
  for( i in 1:nrow(final.configs)){
    
    ### generalized_c_step may not need all of these to be passed in...
    min.info <- generalized_c_step(DATA, weighted.deviations = preproc.DATA$weightedZx, row.weights = preproc.DATA$m, col.weights = preproc.DATA$w, observations_subsample = final.configs[i,])
    final_determinants[i] <- min.info$minimum_determinant
    final_subsamples[i,] <- min.info$observations_subsample
  }
  
  best.order <- order(final_determinants)
  final_determinants <- final_determinants[best.order]
  final_subsamples <- final_subsamples[best.order,]
  
  return( list(final_determinants = final_determinants, final_subsamples = final_subsamples) )
  
}

# formerly cat.c.step

#' @title C step of the GMCD
#' 
#' @description Performs the C ("concentration") step in the search for a sample of observations that produces a minimum determinant from a covariance matrix
#' 
#' @param DATA a data matrix (of presumably all transformed data)
#' @param weighted.deviations the matrix of weighted deviations from \code{\link{ca_preproc}}, specifically \code{$weightedZx}
#' @param row.weights the vector of row weights from \code{\link{ca_preproc}}, specifically \code{$m}
#' @param col.weights the vector of column weights from \code{\link{ca_preproc}}, specifically \code{$w}
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
generalized_c_step <- function(DATA, weighted.deviations, row.weights, col.weights, observations_subsample, max.iters=100, tol=.Machine$double.eps){

  old_determinant <- Inf
  old_center <- NaN
  old_loadings <- matrix(NaN,0,0)
  new_subsample <- old_subsample <- observations_subsample

  for(i in 1:max.iters){
    sub.data <- weighted.deviations[new_subsample,]
    new_center <- colMeans(sub.data)
    svd.res <- tryCatch( {tolerance_svd(sub.data, tol = tol)}, error=function(x) "generalized_c_step: catching error returned from 'tolerance_svd': 'tolerance_svd(sub.data, tol = tol)' failed.") ## this should probably be explained a bit!

    ### should probably use a better catch of thee class here, and should also port this over to continuous_mcd; although that one is less susceptible to these issues
    if(length(svd.res)==3){
      new_determinant <- exp(mean(log(svd.res$d^2)))


      if( (new_determinant <= old_determinant) & (!isTRUE(all.equal(new_determinant,0,tolerance=tol))) ){
        if( center.sigma_checker(old_center, new_center, old_loadings, svd.res$v,tol=tol) & isTRUE(all.equal(new_determinant, old_determinant, tolerance= tol)) ){
          return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
        }
        else{

          old_determinant <- new_determinant
          old_center <- new_center
          old_loadings <- svd.res$v
          old_subsample <- new_subsample

          sup.scores.and.dists <- generalized_scores_dists(DATA, row.weights, col.weights, svd.res$v, svd.res$d)
          new_subsample <- sort(order(sup.scores.and.dists$projected_mahal_dists)[1:length(observations_subsample)]	)

          if( isTRUE(all.equal(sort(new_subsample),sort(old_subsample))) ){
            return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
          }
        }
      }else{
        return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
      }
    }else{
      return( list(observations_subsample = old_subsample, minimum_determinant = old_determinant) )
    }
  }
  ## this is basically the max out which implies the most recent is the best.
  return( list(observations_subsample = new_subsample, minimum_determinant = new_determinant) )
}



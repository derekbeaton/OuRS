## I need to rename some things ("target.data") and also work out the needs of the generalized_ and friends approaches
  #### I need to see existing pipeline code because those may use profiles or weightedZx... I don't recall right now
  #### previous note indicates "## here target.data has to be the weighted deviations (weightedZx)" which makes sense.
    ### I suppose it *could* be Zx too, but the weightedZx is where we get the MD anyways.

#cont.corrmax <- function(population.data,sample.data){

#' @title CorrMax transformation for continuous data
#'
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#'
#' @param target.data a
#' @param center a
#' @param scale a
#' @param loadings a
#' @param singular.values a
#'
#'
#' @author Derek Beaton
#' @export
#'

continuous_corrmax <- function(target.data, center=T, scale=F, loadings, singular.values){

  ## catch NAs
  ## stop if anything is NA; ask that they handle NAs outside of here
  if( any(is.na(target.data)) | any(is.infinite(target.data)) | any(is.nan(target.data)) | any(is.null(target.data)) ){
    stop("continuous_corrmax: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }

  ## should do some error catching here.

  target.data <- ours_scale(target.data, center=center,scale=scale)

  if(missing(loadings) | missing(singular.values)){

    svd.res <- tolerance_svd(target.data)
    loadings <- svd.res$v
    singular.values <- svd.res$d
    svd.res <- NULL ## get rid of it

  }

  res <- corrmax_core_tranform(target.data, loadings, singular.values)
  class(res) <- append(class(res), "continuous")

  res
}

#' @title CorrMax transformation for categorical data
#'
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#'
#' @param target.data a
#' @param loadings a
#' @param singular.values a
#'
#'
#' @author Derek Beaton
#' @export
#'
categorical_corrmax <- function(target.data, loadings, singular.values){

  ## catch NAs
  if( any(is.na(target.data)) | any(is.infinite(target.data)) | any(is.nan(target.data)) | any(is.null(target.data)) ){
    stop("categorical_corrmax: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }

  target.data <- disjunctive_coding(target.data)
  res <- generalized_corrmax(target.data, loadings, singular.values)
  class(res) <- append(class(res), "categorical")

  res
}


#' @title CorrMax transformation for ordinal data
#'
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#'
#' @param target.data a
#' @param mins a
#' @param maxs a
#' @param loadings a
#' @param singular.values a
#'
#'
#' @author Derek Beaton
#' @export
#'
ordinal_corrmax <- function(target.data, mins, maxs, loadings, singular.values){

  ## catch NAs
  if( any(is.na(target.data)) | any(is.infinite(target.data)) | any(is.nan(target.data)) | any(is.null(target.data)) ){
    stop("ordinal_corrmax: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }

  target.data <- thermometer_coding(target.data, mins, maxs)
  res <- generalized_corrmax(target.data, loadings, singular.values)
  class(res) <- append(class(res), "ordinal")

  res
}


#' @title CorrMax transformation for mixed data
#'
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#'
#' @param target.data a
#' @param column.type a
#' @param loadings a
#' @param singular.values a
#'
#'
#' @author Derek Beaton
#' @export
#'
mixed_data_corrmax <- function(target.data, column.type = rep("x", ncol(target.data)), loadings, singular.values){

  ## catch NAs
  if( any(is.na(target.data)) | any(is.infinite(target.data)) | any(is.nan(target.data)) | any(is.null(target.data)) ){
    stop("mixed_data_corrmax: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }

  target.data <- mixed_data_coding(target.data, column.type = column.type)

  res <- generalized_corrmax(target.data, loadings, singular.values)
  class(res) <- append(class(res), "mixed")

  res

}

#' @title CorrMax transformation for correspondence analysis-like data
#'
#' @description A transformation of the data in a row-wise fashion that preserves the Mahalanobis distances yet also provides information on which variables likely contribute to the Mahalanobis distance of the observation
#'
#' @details \code{target.data} are assumed to be correspondence analysis (CA)-like data (i.e., counts). This function performs CA preprocessing before applying the CorrMax procedure.
#' This particular function is the primary function for \code{\link{categorical_corrmax}}, \code{\link{ordinal_corrmax}}, and \code{\link{mixed_data_corrmax}}.
#'
#' @param target.data a
#' @param loadings a
#' @param singular.values a
#'
#'
#' @author Derek Beaton
#' @export
#'
generalized_corrmax <- function(target.data, loadings, singular.values){

  ## should do some error catching here.

  ## catch NAs
  if( any(is.na(target.data)) | any(is.infinite(target.data)) | any(is.nan(target.data)) | any(is.null(target.data)) ){
    stop("generalized_corrmax: NA, Inf, -Inf, NULL, and NaN are not allowed.")
  }

  preproc.DATA <- ca_preproc(target.data, compact = T) ## this could be more efficient...
  if(missing(loadings) | missing(singular.values)){

    svd.res <- tolerance_svd(preproc.DATA$weightedZx)
    loadings <- svd.res$v
    singular.values <- svd.res$d
    svd.res <- NULL ## get rid of it

  }
  res <- corrmax_core_tranform(preproc.DATA$weightedZx, loadings, singular.values)
  class(res) <- append(class(res), "generalized")

  res
}





#' @noRd
#'
#' @title The core of the CorrMax functions
#'
#' @description something
#'
#' @param target.data a
#' @param loadings a
#' @param singular.values a
#'
#' @author Derek Beaton

### I guess technically this is a "covmax"
corrmax_core_tranform <- function(target.data, loadings, singular.values){

  ## should do some error catching here.


  diag.sampcov <- 1/sqrt(diag(crossprod(t(loadings) * singular.values)))
  inv.DSD.half <- GSVD::invsqrt_psd_matrix( tcrossprod(t(t(loadings) * singular.values) * diag.sampcov) )
  transformed_data <- target.data %*% (t(inv.DSD.half) * diag.sampcov)

  contributions <- transformed_data^2
  mahal_dists <- rowSums(contributions)

  res <- list(
    diag_sample_cov = diag.sampcov,
    corr_max_transform = inv.DSD.half,
    transformed_data = transformed_data,
    contributions = contributions,
    mahal_dists = mahal_dists,
    percentage_contributions = sweep(contributions,1,mahal_dists,"/")*100
  )

  class(res) <- c("list", "OuRS", "CorrMax")
  res
}

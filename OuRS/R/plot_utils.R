### plot utils

### the component plot is not correctly scaled according to D, but, that's OK this is just a simple visual
#' @title Plot components from PCA or CA
#' @description a simple 2D plot for observation (or variable) component scores or vectors to visualize PCA or CA
#' @param scores a numeric matrix with the scores (fi or fj) or vectors (u/p or v/q)
#' @param axes a numeric vector of length 2. Indicates which axes (columns) of \code{scores} to plot
#' @param pch see \code{\link{plot}} and \code{\link{par}}. Default here is '20'
#' @param col see \code{\link{plot}} and \code{\link{par}}. Default here is 'mediumorchid4'
#' @param main see \code{\link{plot}} and \code{\link{par}}. Default here is 'Component scores'
#' @param xlab see \code{\link{plot}} and \code{\link{par}}. Default here is '\code{paste0("Component ",axes[1])}'
#' @param ylab see \code{\link{plot}} and \code{\link{par}}. Default here is '\code{paste0("Component ",axes[2])}'
#' @param xlim see \code{\link{plot}} and \code{\link{par}}. Default here is '\code{c(-max(abs(scores[,axes])),max(abs(scores[,axes])))*1.3}'
#' @param ylim see \code{\link{plot}} and \code{\link{par}}. Default here is '\code{c(-max(abs(scores[,axes])),max(abs(scores[,axes])))*1.3}'
#' @param asp see \code{\link{plot}} and \code{\link{par}}. Default here is '1'
#' @param pos see \code{\link{plot}} and \code{\link{par}}. Default here is '3'
#' @param display_names a logical (boolean). If \code{TRUE} plot the row names of \code{scores}
#' @param cex see \code{\link{plot}} and \code{\link{par}}. Default here is '1'
#' @param text.cex like \code{cex} but for text plotting when \code{display_names = TRUE}. See also \code{\link{plot}} and \code{\link{par}}. Default here is '1'
#' @param ... Additional arguments to be passed for plotting. See \code{\link{plot}} and \code{\link{par}}.
#' 
#' @author Derek Beaton
#' @export
component_plot <- function(scores, axes=c(1,2), pch=20, col="mediumorchid4",
                           main="Component scores",
                           xlab=paste0("Component ",axes[1]),
                           ylab=paste0("Component ",axes[2]),
                           xlim=c(-max(abs(scores[,axes])),max(abs(scores[,axes])))*1.3,
                           ylim=c(-max(abs(scores[,axes])),max(abs(scores[,axes])))*1.3,
                           asp=1, pos=3, display_names=T,cex=1,text.cex=1,
                           ...){

  plot(0, type="n", xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, axes=F, asp=asp)
  abline(h=0,lty=2,lwd=2, col="grey60")
  abline(v=0,lty=2,lwd=2, col="grey60")
  points(scores[,axes], col=col, pch=pch, cex=cex, ...)

  if (display_names) {
    text(scores[, axes], labels = rownames(scores),
         pos = pos, col = col, cex = text.cex)
  }

}


## pass in where we want the cutoffs
#' @title Distance-distance plot
#' @description A plot of standard vs. robust Mahalanobis distances (with optional transformations)
#' @param ours_mcd_list an OuRS MCD class object (of type list). Returned from \code{continuous_mcd}, \code{categorical_mcd}, \code{ordinal_mcd}, or \code{mixed_data_mcd}
#' @param md_cutoff numeric. A value for Mahalanobis distances to display a cutoff. If invalid, \code{quantile(x, probs=.95)} used
#' @param robust_md_cutoff numeric. A value for robust Mahalanobis distances to display a cutoff. If invalid, \code{quantile(x, probs=.95)} used
#' @param dist_transform character of type "none", "sqrt", or "log". If "sqrt" the distances will be transformed by square root, if "log" the distances will be transformed by the natural log. If "none" , no transformation will be performed.

dd_plot <- function(ours_mcd_list, md_cutoff = NA, robust_md_cutoff = NA,  dist_transform = "none"){
  
  if(!inherits(ours_mcd_list,c("list", "OuRS", "MCD"))){
    stop("dd_plot: 'ours_mcd_list' is not a recognized OuRS MCD object.")
  }
  
  ## 3 options for transforms
  if(length(dist_transform) != 1){
    stop("dd_plot: 'dist_transform' must be of length 1")
  }
  if( !(dist_transform %in% c("none","sqrt","log")) ){
    stop("dd_plot: 'dist_transform' is not one of the recognized types ('none','sqrt','log')")
  }
  
  
  # if a transform for Ds, also transform h & v
  xy <- cbind(ours_mcd_list$dists$mahal_dists, ours_mcd_list$dists$robust_mahal_dists)
  
  # if(!is.numeric(md_cutoff) | is.na(md_cutoff) | is.nan(md_cutoff) | is.infinite(md_cutoff) | is.null(md_cutoff)){
  #   md_cutoff <- quantile(xy[,1], probs = .95)
  # }
  # if(md_cutoff < min(xy) | md_cutoff > max(xy)){
  #   md_cutoff <- quantile(xy[,1], probs = .95)
  # }
  # 
  # if(!is.numeric(robust_md_cutoff) | is.na(robust_md_cutoff) | is.nan(robust_md_cutoff) | is.infinite(robust_md_cutoff) | is.null(robust_md_cutoff)){
  #   robust_md_cutoff <- quantile(xy[,1], probs = .95)
  # }
  # if(robust_md_cutoff < min(xy) | robust_md_cutoff > max(xy)){
  #   robust_md_cutoff <- quantile(xy[,1], probs = .95)
  # }
  
  
  
  if(dist_transform=="sqrt"){
    
    plot(sqrt(xy), xlab = "Square root of squared Mahalanobis Distances", ylab = "Square root of squared robust Mahalanobis Distances", main = "distance-distance plot", pch = 21, bg = "grey80", col="black")
    if( is.numeric(robust_md_cutoff )){
      if( !is.na(robust_md_cutoff) & !is.infinite(robust_md_cutoff) & is.nan(robust_md_cutoff) ){
        abline(h = sqrt(robust_md_cutoff), col="firebrick3")
      }
    }
    if( is.numeric(md_cutoff )){
      if( !is.na(md_cutoff) & !is.infinite(md_cutoff) & is.nan(md_cutoff) ){
        abline(v = sqrt(md_cutoff), col="steelblue4")
      }
    }
    
  }else if(dist_transform=="log"){
    
    plot(log(xy), xlab = "Natural log of squared Mahalanobis Distances", ylab = "Natural log of squared robust Mahalanobis Distances", main = "distance-distance plot", pch = 21, bg = "grey80", col="black")
    if( is.numeric(robust_md_cutoff )){
      if( !is.na(robust_md_cutoff) & !is.infinite(robust_md_cutoff) & is.nan(robust_md_cutoff) ){
        abline(h = log(robust_md_cutoff), col="firebrick3")
      }
    }
    if( is.numeric(md_cutoff )){
      if( !is.na(md_cutoff) & !is.infinite(md_cutoff) & is.nan(md_cutoff) ){
        abline(v = log(md_cutoff), col="steelblue4")
      }
    }
    
  }else{
    
    plot(xy, xlab = "Squared Mahalanobis Distances", ylab = "Squared robust Mahalanobis Distances", main = "distance-distance plot", pch = 21, bg = "grey80", col="black")
    if( is.numeric(robust_md_cutoff )){
      if( !is.na(robust_md_cutoff) & !is.infinite(robust_md_cutoff) & is.nan(robust_md_cutoff) ){
        abline(h = robust_md_cutoff, col="firebrick3")
      }
    }
    if( is.numeric(md_cutoff )){
      if( !is.na(md_cutoff) & !is.infinite(md_cutoff) & is.nan(md_cutoff) ){
        abline(v = md_cutoff, col="steelblue4")
      }
    }
    
  }
  
  
}


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
#' @param ours_mcd_list
#' @param horizontal_line
#' @param vertical_cutoff
#' @param dist_transform
#' 
dd_plot <- function(ours_mcd_list, horizontal_line = NA, vertical_cutoff = NA, dist_transform = "none"){
  
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
  if(dist_transform=="sqrt"){
    
    plot(sqrt(xy), )
    
  }else if(dist_transform=="log"){
    
    plot(log(xy), )
    
  }else{
    
    plot(xy, )
    
  }
  
  
}


#' @title Make tolerance ellipse for a 2D or x-y scatter plot
#' @description Compute the tolerance ellipse for x-y data
#' @details Note: This function is effectively a wrapper for a set of functions copied from the R package SIBER.
#' Here, many of the SIBER functions are private to \code{tolerance_ellipse}, but were copied directly from SIBER 2.1.0. 
#' So users of \code{OuRS} only interface with \code{tolerance_ellipse} and not any of the private functions
#'
#' @param DATA a numeric matrix with (presumably) two columns. This function will only make use of the first two.
#' @param ellipse.alpha numeric in the range of (.5,1). The percentage for a tolerance ellipse
#' @param mcd.alpha numeric in the range of (.5,1).  The percentage for a robust tolerance ellipse (by way of the MCD)
#' @param xlab see \code{\link{plot}} and \code{\link{par}}. Default here is '\code{colnames(DATA)[1]}'
#' @param ylab see \code{\link{plot}} and \code{\link{par}}. Default here is '\code{colnames(DATA)[2]}'
#' @param graphs logical (boolean). Default is \code{FALSE}. When \code{FALSE} ellipses are added to an existing plot, when \code{TRUE} a new x-y plot is made with ellipses
#' 
#' @author Derek Beaton for \code{tolerance_ellipse}; Andrew Jackson and Andrew Parnell for all functions inside (i.e., \code{pointsToEllipsoid}, \code{ellipsoidTransform}, \code{ellipseInOut}, \code{addEllipse}, \code{genCircle})
#' @seealso https://CRAN.R-project.org/package=SIBER
#' @export
tolerance_ellipse <- function(DATA, ellipse.alpha=.75, mcd.alpha=.75, xlab=colnames(DATA)[1], ylab=colnames(DATA)[2], graphs=F){

  ## private function. STOLEN FROM SIBER 2.1.0
  pointsToEllipsoid <- function (X, Sigma, mu)
  {
    if (ncol(Sigma) != nrow(Sigma))
      stop("Sigma must be a square matrix")
    if (ncol(X) != ncol(Sigma))
      stop("number of columns in X must \n                                  be of same dimension as Sigma")
    if (length(mu) != ncol(Sigma))
      stop("length of mu must \n                                  be of same dimension as Sigma")
    eig <- eigen(Sigma)
    SigSqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
    Z <- t(apply(X, 1, ellipsoidTransform, SigSqrt, mu))
    return(Z)
  }

  ## pTE private function. STOLEN FROM SIBER 2.1.0
  ellipsoidTransform <- function (x, SigSqrt, mu)
  {
    return(solve(SigSqrt, x - mu))
  }

  ## private function. STOLEN FROM SIBER 2.1.0
  ellipseInOut <- function (Z, p = 0.95, r = NULL)
  {
    if (is.null(r)) {
      r <- stats::qchisq(p, df = ncol(Z))
    }
    inside <- rowSums(Z^2) < r
    return(inside)
  }

  ## private function. STOLEN FROM SIBER 2.1.0
  addEllipse <- function (mu, sigma, m = NULL, n = 100, p.interval = NULL, ci.mean = FALSE,small.sample = FALSE, do.plot = TRUE, ...)
  {
    if (small.sample & is.null(m))
      message("A sample size number given by m is \n required when small.sample is TRUE")
    if (ci.mean & is.null(m))
      message("A sample size number given by m is \n  required when plotting confidence \n ellipses of the mean with ci.mean is TRUE")
    ifelse(ci.mean, c.scale <- m, c.scale <- 1)
    ifelse(small.sample, q <- (m - 1)/(m - 2), q <- 1)
    ifelse(is.null(p.interval), r <- 1, r <- sqrt(stats::qchisq(p.interval,
                                                                df = 2)))
    e = eigen(sigma/c.scale)
    SigSqrt = e$vectors %*% diag(sqrt(e$values * q)) %*% t(e$vectors)
    cc <- genCircle(n, r)
    back.trans <- function(x) {
      return(SigSqrt %*% x + mu)
    }
    ML.ellipse = t(apply(cc, 1, back.trans))
    if (grDevices::dev.cur() > 1 & do.plot) {
      graphics::lines(ML.ellipse, ...)
    }
    return(ML.ellipse)
  }

  ## private function. STOLEN FROM SIBER 2.1.0
  genCircle <- function (n = 1000, r)
  {
    theta = seq(0, 2 * pi, length = n)
    x = r * cos(theta)
    y = r * sin(theta)
    return(cbind(x, y))
  }

  if(alpha < .5){
    alpha <- .5
  }
  if(alpha > 1){
    alpha <- 1
  }
  
  mcd <- covMcd(DATA,alpha = mcd.alpha)
  mcd.center <- mcd$center
  mcd.cov <- mcd$cov
  data.center <- colMeans(DATA)
  data.cov <- cov(DATA)



  if(graphs){
    x1 <- c(-max(abs(DATA[,1]))*.05,max(abs(DATA[,1])))*1.1
    y1 <- c(-max(abs(DATA[,2]))*.05,max(abs(DATA[,2])))*1.1


    plot(DATA, xlim = x1, ylim = y1,pch=20,col="grey80", main="", xlab=xlab, ylab=ylab,cex=.5)

    rob.ellipse <- addEllipse(mcd.center,mcd.cov,p.interval = ellipse.alpha,col="blue",lty=2)
    classic.ellipse <- addEllipse(data.center,data.cov,p.interval = ellipse.alpha,col="red",lty=2)

    abline(v=max(rob.ellipse[,1]), h=max(rob.ellipse[,2]),lty=1,col="blue")
    points(DATA[which(DATA[,1] >= max(rob.ellipse[,1]) | DATA[,2] >= max(rob.ellipse[,2])),],bg="blue",pch=21,cex=1)
    abline(v=max(classic.ellipse[,1]), h=max(classic.ellipse[,2]),lty=1,col="red")
    points(DATA[which(DATA[,1] >= max(classic.ellipse[,1]) | DATA[,2] >= max(classic.ellipse[,2])),],bg="red",pch=21,cex=2)
    legend("bottomright",legend=c(paste0("Classic ellipse alpha = ", ellipse.alpha),paste0("Robust ellipse alpha = ", ellipse.alpha, "with MCD alpha = ", mcd.alpha)), col=c("red","blue"), lty=c(2,1))
  }else{
    rob.ellipse <- addEllipse(mcd.center,mcd.cov,p.interval = ellipse.alpha,col="blue",lty=2,do.plot = F)
    classic.ellipse <- addEllipse(data.center,data.cov,p.interval = ellipse.alpha,col="red",lty=2,do.plot = F)
  }

  return(
    list(
      x.robust.cutoff=max(rob.ellipse[,1]),
      x.classic.cutoff=max(classic.ellipse[,1]),
      y.robust.cutoff=max(rob.ellipse[,2]),
      y.classic.cutoff=max(classic.ellipse[,2]),

      x.robust.outliers = (DATA[,1] >= max(rob.ellipse[,1])),
      x.classic.outliers= (DATA[,1] >= max(classic.ellipse[,1])),
      y.robust.outliers = (DATA[,2] >= max(rob.ellipse[,2])),
      y.classic.outliers= (DATA[,2] >= max(classic.ellipse[,2]))
    )
  )

}


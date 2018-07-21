tol.ellipse <- function(dat,ellipse.alpha=.75,mcd.alpha=.75,xlab=colnames(dat)[1],ylab=colnames(dat)[2],graphs=F){

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

  mcd <- covMcd(dat,alpha = mcd.alpha)
  mcd.center <- mcd$center
  mcd.cov <- mcd$cov
  data.center <- colMeans(dat)
  data.cov <- cov(dat)



  if(graphs){
    x1 <- c(-max(abs(dat[,1]))*.05,max(abs(dat[,1])))*1.1
    y1 <- c(-max(abs(dat[,2]))*.05,max(abs(dat[,2])))*1.1


    plot(dat, xlim = x1, ylim = y1,pch=20,col="grey80", main="", xlab=xlab, ylab=ylab,cex=.5)

    rob.ellipse <- addEllipse(mcd.center,mcd.cov,p.interval = ellipse.alpha,col="blue",lty=2)
    classic.ellipse <- addEllipse(data.center,data.cov,p.interval = ellipse.alpha,col="red",lty=2)

    abline(v=max(rob.ellipse[,1]), h=max(rob.ellipse[,2]),lty=1,col="blue")
    points(dat[which(dat[,1] >= max(rob.ellipse[,1]) | dat[,2] >= max(rob.ellipse[,2])),],bg="blue",pch=21,cex=1)
    abline(v=max(classic.ellipse[,1]), h=max(classic.ellipse[,2]),lty=1,col="red")
    points(dat[which(dat[,1] >= max(classic.ellipse[,1]) | dat[,2] >= max(classic.ellipse[,2])),],bg="red",pch=21,cex=2)
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

      x.robust.outliers = (dat[,1] >= max(rob.ellipse[,1])),
      x.classic.outliers= (dat[,1] >= max(classic.ellipse[,1])),
      y.robust.outliers = (dat[,2] >= max(rob.ellipse[,2])),
      y.classic.outliers= (dat[,2] >= max(classic.ellipse[,2]))
    )
  )

}

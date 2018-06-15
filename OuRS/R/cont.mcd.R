	## I do not know what robust.mahals is as a flag...
#cont.mcd <- function(data, center=T, scale=F, collinearity.stop=T, alpha=.75, robust.mahals=T, h.size.abs=F, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){
cont.mcd <- function(data, center=T, scale=F, collinearity.stop=T, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){


  if(ncol(data) > (nrow(data)*.9)){
  	stop("Column:Row ratio is too high. Please use another method (e.g., robPCA, rPCA, POET).")
  }

  ## collinearity check
  svd.res <- svd(expo.scale(data,center=center,scale=scale))
  collinear.components <- which(svd.res$d^2 < tol)
  if(length(collinear.components)>0){
    if(collinearity.stop){
      stop("Data are collinear. Some variables are likely combinations of other variables. Please perform a plain SVD and inspect the small value singular values and their respective vectors.")
    }else{
      warning("Data are collinear. Some variables are likely combinations of other variables. You have chosen to ignore these components for the MCD computation. MCD will proceed.")
    }
  }


  ## sample finder
  mcd.samples <- cont.mcd.find.sample(data, center=center, scale=scale, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent,tol=tol)
  ## only grab the top sample.
  best.sample <- mcd.samples$final.orders[1,]

  ### NOTE: I should pass the tol parameter through to all places where tolerance.svd gets called.
  ## compute standard distances (especially Mahal)
  tsvd.res <- tolerance.svd(expo.scale(data,center=center,scale=scale),tol = tol)
  mahals <- rowSums(tsvd.res$u^2)
  chis <- rowSums( (tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T))^2 )

  ## get robust mean & cov (loadings)
  rob.sample <- expo.scale(data[best.sample,],center=center,scale=scale)
  rob.center <- attributes(rob.sample)$`scaled:center`
  rob.scale <- attributes(rob.sample)$`scaled:scale`
  robust.tsvd.res <- tolerance.svd(rob.sample,tol = tol)

  ## call to function that computes robust mahal
  robust.dists <- cont.sup.fi.u(data,center=rob.center,scale=rob.scale,robust.tsvd.res$v,robust.tsvd.res$d)
  robust.mahals <- rowSums(robust.dists$sup.u^2)
  robust.chis <- rowSums(robust.dists$sup.fi^2)


  # compute ODs.
  if(ncol(tsvd.res$u)==ncol(robust.dists$sup.u)){
    u.od <- rowSums((tsvd.res$u - robust.dists$sup.u)^2)
  }else{
    u.od <- NA #for now
  }

  if(ncol(tsvd.res$u)==ncol(robust.dists$sup.fi)){
    fi.od <- rowSums(( sweep(tsvd.res$u,2,tsvd.res$d,"*") - robust.dists$sup.fi)^2)  # should be identical...
  }else{
    fi.od <- NA #for now
  }

  # return(
  #   list( best.det=mcd.samples$final.dets[1],
  #         best.sample=best.sample,
  #         best.loadings= robust.tsvd.res$v,
  #         best.svs=robust.tsvd.res$d,
  #         best.rob.md = robust.mahals,
  #         best.rob.chid= robust.chis,
  #         best.m.od= mahal.od,
  #         best.chi.od= chi.od,
  #         md=mahals,
  #         chid=chis,
  #         final.sets = list(final.dets = mcd.samples$final.dets, final.orders = mcd.samples$final.orders))
  # )
  res <- list(
    cov = list(loadings = robust.tsvd.res$v,
               singular.values = robust.tsvd.res$d,
               center = rob.center,
               scale = rob.scale
    ),
    dists = list(rob.md = robust.mahals,
                 rob.chid = robust.chis,
                 md = mahals,
                 chid = chis,
                 u.od = u.od,
                 fi.od = fi.od),
    det.samps = list(dets = mcd.samples$final.dets,
                     samples = mcd.samples$final.orders)
  )
  class(res) <- c("contMCD","list")
  return(res)
}

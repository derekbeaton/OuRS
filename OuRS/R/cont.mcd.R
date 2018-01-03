	## I do not know what robust.mahals is as a flag...
#cont.mcd <- function(data, center=T, scale=F, collinearity.stop=T, alpha=.75, robust.mahals=T, h.size.abs=F, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){
cont.mcd <- function(data, center=T, scale=F, collinearity.stop=T, alpha=.75, h.size.abs=F, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){


  if(ncol(data) > (nrow(data)*.9)){
  	stop("Column:Row ratio is too high. Please use another method (e.g., robPCA, rPCA, POET).")
  }

  ## collinearity check
  svd.res <- svd(expo.scale(data,center=center,scale=scale))
  collinear.components <- which(svd.res$d < tol)
  if(length(collinear.components)>0){
    if(collinearity.stop){
      stop("Data are collinear. Some variables are likely combinations of other variables. Please perform a plain SVD and inspect the small value singular values and their respective vectors.")
    }else{
      warning("Data are collinear. Some variables are likely combinations of other variables. You have chosen to ignore these components for the MCD computation. MCD will proceed.")
    }
  }


  ## sample finder
  mcd.samples <- cont.mcd.find.sample(data, center=center, scale=scale, alpha=alpha, h.size.abs=h.size.abs, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent)
  ## only grab the top sample.
  best.sample <- mcd.samples$final.orders[1,]


  ## compute standard distances (especially Mahal)
  tsvd.res <- tolerance.svd(expo.scale(data,center=center,scale=scale))
  mahals <- rowSums(tsvd.res$u^2)
  chis <- rowSums( (tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T))^2 )

  ## get robust mean & cov (loadings)
  rob.sample <- expo.scale(data[best.sample,],center=center,scale=scale)
  rob.center <- attributes(rob.sample)$`scaled:center`
  rob.scale <- attributes(rob.sample)$`scaled:scale`
  robust.tsvd.res <- tolerance.svd(rob.sample)

  ## call to function that computes robust mahal
  robust.dists <- sup.fi.u(data,center=rob.center,scale=rob.scale,robust.tsvd.res$v,robust.tsvd.res$d)
  robust.mahals <- rowSums(robust.dists$sup.u^2)
  robust.chis <- rowSums(robust.dists$sup.fi^2)


	# compute ODs.
	mahal.od <- rowSums((tsvd.res$u - robust.dists$sup.u)^2)
	chi.od <- rowSums(((tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T)) - robust.dists$sup.fi)^2)

	#
	return(
	  list( best.det=mcd.samples$final.dets[1],
	        best.order=mcd.samples$final.orders[1,],
	        best.loadings= robust.tsvd.res$v,
	        best.center= rob.center,
	        best.scale= rob.scale,
	        best.rob.md = robust.mahals,
	        best.rob.chid= robust.chis,
	        best.m.od= mahal.od,
	        best.chi.od= chi.od,
	        md= mahals,
	        chid=(tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T)),
	        final.sets = list(final.dets = mcd.samples$final.dets, final.orders = mcd.samples$final.orders)) )
}

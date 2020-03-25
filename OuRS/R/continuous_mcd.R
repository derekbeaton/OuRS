
##
#'
#' @export
#'

continuous_mcd <- function(data, center=T, scale=F, collinearity.stop=T, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){


  if(ncol(data) > (nrow(data)*.9)){
  	stop("continuous_mcd: the column-to-row ratio is too high so 'mcd' cannot be performed")
  }

  ## collinearity check
  svd.res <- tolerance_svd(ours_scale(data,center=center,scale=scale), tol = NULL)
  if( any(svd.res$d^2 < tol) ) {
    if(collinearity.stop){
      stop("continuous_mcd: Matrix is singular. It is likely that data are collinear where some variables are combinations of other variables.")
    }
  }

  ## check alpha here.


  ## sample finder
  mcd.samples <- cont_mcd_find_sample(data, center=center, scale=scale, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent,tol=tol)

  ## get sample distances
  tsvd.res <- tolerance_svd(ours_scale(data,center=center,scale=scale),tol = tol)
  mahals <- rowSums(tsvd.res$u^2)
  chis <- rowSums( (tsvd.res$u  * matrix(tsvd.res$d,nrow(tsvd.res$u),ncol(tsvd.res$u),byrow=T))^2 )

  ## get robust mean & cov (loadings)
  rob.sample <- ours_scale(data[mcd.samples$final.orders[1,],],center=center,scale=scale)
  rob.center <- attributes(rob.sample)$`scaled:center`
  rob.scale <- attributes(rob.sample)$`scaled:scale`
  robust.tsvd.res <- tolerance_svd(rob.sample,tol = tol)

  ## get robust distances
  robust.vectors.and.scores <- cont.sup.fi.u(data,center=rob.center,scale=rob.scale,robust.tsvd.res$v,robust.tsvd.res$d)
  robust.mahals <- rowSums(robust.vectors.and.scores$sup.u^2)
  robust.chis <- rowSums(robust.vectors.and.scores$sup.fi^2)


  res <- list(
    cov = list(loadings = robust.tsvd.res$v,
               singular.values = robust.tsvd.res$d,
               center = rob.center,
               scale = rob.scale
    ),
    dists = list(rob.md = robust.mahals,
                 rob.chid = robust.chis,
                 md = mahals,
                 chid = chis),
    det.samps = list(dets = mcd.samples$final.dets,
                     samples = mcd.samples$final.orders)
  )
  class(res) <- c("list", "contMCD")


  res
}

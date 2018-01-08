cat.mcd <- function(data, make.data.disjunctive=F, alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){

  if(make.data.disjunctive){
    data <- make.data.nominal(data)
    make.data.disjunctive <- F
  }

  if(ncol(data) > (nrow(data)*.9)){
    stop("Column:Row ratio is too high. Please use another method (e.g., robPCA, rPCA, POET).")
  }


  ## sample finder
  mcd.samples <- cat.mcd.find.sample(data, make.data.disjunctive=make.data.disjunctive, alpha=alpha, num.subsets=num.subsets, max.total.iters=max.total.iters, top.sets.percent=top.sets.percent)
  ## only grab the top sample.
  best.sample <- mcd.samples$final.orders[1,]

  preproc.data <- ca.preproc(data) ## this could be more efficient...
  #profiles <- (diag(1/preproc.data$m) %*% preproc.data$Ox)
  profiles <- sweep(preproc.data$Ox,1,preproc.data$m,"/")

  ca.res <- ca(data)
  mahals <- rowSums(ca.res$u^2)
  chis <- rowSums(ca.res$fi^2)


  ## get robust mean & cov (loadings)
  rob.sample <- preproc.data$weightedZx[best.sample,]
  robust.tsvd.res <- tolerance.svd(rob.sample)

  ## call to function that computes robust mahal
  robust.dists <- cat.sup.fi.u(profiles, preproc.data$m, preproc.data$w, robust.tsvd.res$v, robust.tsvd.res$d)
  robust.mahals <- rowSums(robust.dists$sup.u^2)
  robust.chis <- rowSums(robust.dists$sup.fi^2)

  # compute ODs.
  if(ncol(ca.res$u)==ncol(robust.dists$sup.u)){
    u.od <- rowSums((ca.res$u - robust.dists$sup.u)^2)
  }else{
    u.od <- NA #for now
  }

  if(ncol(ca.res$fi)==ncol(robust.dists$sup.fi)){
    fi.od <- rowSums((ca.res$fi - robust.dists$sup.fi)^2)  # should be identical...
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
  return(list(
  cov = list(loadings = robust.tsvd.res$v,
              singular.values = robust.tsvd.res$d
  ),
  dists = list(rob.md = robust.mahals,
               rob.chid = robust.chis,
               md = mahals,
               chid = chis,
               u.od = u.od,
               fi.od = fi.od),
  det.samps = list(dets = mcd.samples$final.dets,
                   samples = mcd.samples$final.orders)
  ))

}

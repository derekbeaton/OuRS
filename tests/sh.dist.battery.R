
sh.dist.battery <- function(DATA,sh.pca.res,center=T,scale=F,comps=NULL){

  full.comps <- F
  if(is.null(comps)){
    comps <- 1:max(dim(sh.pca.res$loadings.cors)[1:2])
    full.comps <- T
  }

  if(length(setdiff(1:max(dim(sh.pca.res$loadings.cors)[1:2]),comps))==0){
    full.comps <- T
  }

  comps <- sort(comps)

  #pred.fi.mean <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,mean)
  pred.fi.median <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,median)
  pred.fi.max <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,max)
  pred.fi.min <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,min)
  pred.fi.iqr <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,IQR)
  #pred.fi.sd <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,sd)
  #pred.fi.diff <- apply(apply(sh.pca.res$pred.fi.array^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})

  #pred.u.mean <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,mean)
  pred.u.median <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,median)
  pred.u.max <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,max)
  pred.u.min <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,min)
  pred.u.iqr <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,IQR)
  #pred.u.sd <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,sd)
  #pred.u.diff <- apply(apply(sh.pca.res$pred.u.array^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})

  min.shr.det <- which(sh.pca.res$sh.dets == min(sh.pca.res$sh.dets), arr.ind = T)
  if (min.shr.det[2] == 1) {
    shr.min.det.sample <- sh.pca.res$sh1.orders[min.shr.det[1],]
  } else {
    shr.min.det.sample <- sh.pca.res$sh2.orders[min.shr.det[1],]
  }

  shr.c.step <- cont.c.step(data = DATA, obs.order = shr.min.det.sample, center = center, scale = scale, max.iters = 100)
  shr.c.step.pca <- tolerance.svd(scale(DATA[shr.c.step$obs.order, ], center = center, scale = scale))

  scaled.dat <- expo.scale(DATA[shr.c.step$obs.order, ], center = center, scale = scale)
  shr.c.step.pred.fi <- scale(DATA, center = attributes(scaled.dat)$`scaled:center`, scale = attributes(scaled.dat)$`scaled:scale`) %*% shr.c.step.pca$v
  shr.c.step.pred.u <- scale(DATA, center = attributes(scaled.dat)$`scaled:center`, scale = attributes(scaled.dat)$`scaled:scale`) %*% shr.c.step.pca$v %*% diag(1/shr.c.step.pca$d)
  shr.c.step.score <- rowSums(shr.c.step.pred.fi^2)
  shr.c.step.mahal <- rowSums(shr.c.step.pred.u^2)


  if(!full.comps){
    #lr.pred.fi.mean <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,mean)
    lr.pred.fi.median <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,median)
    lr.pred.fi.max <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,max)
    lr.pred.fi.min <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,min)
    lr.pred.fi.iqr <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,IQR)
    #lr.pred.fi.sd <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,sd)
    #lr.pred.fi.diff <- apply(apply(sh.pca.res$pred.fi.array[,comps,]^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})

    #lr.pred.u.mean <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,mean)
    lr.pred.u.median <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,median)
    lr.pred.u.max <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,max)
    lr.pred.u.min <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,min)
    lr.pred.u.iqr <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,IQR)
    #lr.pred.u.sd <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,sd)
    #lr.pred.u.diff <- apply(apply(sh.pca.res$pred.u.array[,comps,]^2,c(1,3),sum),1,function(x){abs(max(x)-min(x))})

    lr.scaled.dat <- expo.scale(DATA, center = center, scale = scale)
    svd.res <- tolerance.svd(lr.scaled.dat)

    low.rank.fi <- svd.res$u[, comps] %*% diag(svd.res$d[comps])
    low.rank.DATA <- low.rank.fi %*% t(svd.res$v[, comps])
    low.rank.DATA <- low.rank.DATA * matrix(attributes(lr.scaled.dat)$`scaled:scale`, nrow(low.rank.DATA), ncol(low.rank.DATA), byrow = T) + matrix(attributes(lr.scaled.dat)$`scaled:center`, nrow(low.rank.DATA), ncol(low.rank.DATA), byrow = T)
    od <- rowSums((DATA - low.rank.DATA)^2)


    low.rank.sd <- rowSums(low.rank.fi^2)
    low.rank.md <- rowSums(svd.res$u[, comps]^2)


    lr.shr.c.step <- cont.c.step(data = low.rank.DATA, obs.order = shr.min.det.sample, center = center, scale = scale, max.iters = 100)
    lr.shr.c.step.pca <- tolerance.svd(scale(low.rank.DATA[lr.shr.c.step$obs.order, ], center = center, scale = scale))
    # alt od
    # T = (DATA - ROB.center) * ROB.EIG.VECS
    # (DATA - ROB.center) - (T * ROB.EIG.VECS')
        ## so (T * ROB.EIG.VECS) is some sort of rebuilt robust data...


    lr.scaled.sub.dat <- expo.scale(low.rank.DATA[lr.shr.c.step$obs.order, ],center = center, scale = scale)

    lr.shr.c.step.pred.fi <- scale(DATA, center = attributes(lr.scaled.sub.dat)$`scaled:center`, scale = attributes(lr.scaled.sub.dat)$`scaled:scale`) %*% lr.shr.c.step.pca$v
    lr.shr.c.step.pred.u <- scale(DATA, center = attributes(lr.scaled.sub.dat)$`scaled:center`, scale = attributes(lr.scaled.sub.dat)$`scaled:scale`) %*% lr.shr.c.step.pca$v %*% diag(1/lr.shr.c.step.pca$d)
    lr.shr.c.step.score <- rowSums(lr.shr.c.step.pred.fi^2)
    lr.shr.c.step.mahal <- rowSums(lr.shr.c.step.pred.u^2)


    return(list(
                 dists=cbind(pred.fi.median, pred.fi.iqr, pred.fi.max, pred.fi.min,
                             pred.u.median,pred.u.iqr, pred.u.max, pred.u.min,
                              shr.c.step.score,shr.c.step.mahal),

                 lr.sh.dists=cbind(lr.pred.fi.median, lr.pred.fi.iqr, lr.pred.fi.max, lr.pred.fi.min,
                                   lr.pred.u.median, lr.pred.u.iqr, lr.pred.u.max, lr.pred.u.min,
                                   lr.shr.c.step.score,lr.shr.c.step.mahal),

                 lr.dists=cbind(low.rank.sd,low.rank.md,od)

                 ))

  }

  return(list(dists=cbind(pred.fi.median,pred.fi.iqr, pred.fi.max, pred.fi.min,
                          pred.u.median, pred.u.iqr, pred.u.max, pred.u.min,
                          shr.c.step.score,shr.c.step.mahal)))

}

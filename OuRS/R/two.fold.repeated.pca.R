two.fold.repeated.pca <- function(DATA,center=T,scale=F,iters=500,sh1.size=.5,k=0){


  if(sh1.size > .9 | sh1.size < .1){
    warning("sh1.size is greater than or equal to 90% or less than or equal to 10%. Setting sh1.size to .5")
    sh1.size <- .5
  }

  if(nrow(DATA) < 20){
    warning("DATA has less than 20 rows. sh1.size will automatically be set to 50%")
    sh1.size <- .5
  }

  ## do the intitial PCA here anyways. We'll be able to use it for a lot of these things.
  #pca.res <- gsvd(expo.scale(DATA,center=center,scale=scale),k=k)
  pca.res <- pca(DATA,center=center,scale=scale,k=k,compact=F)
  max.rank <- length(pca.res$d.orig)
    ## I'm not really using the original PCA much now but can later for lots of things, e.g., comparing predictions against this.

  sh1.orders <- matrix(NA,iters,ceiling(nrow(DATA)*sh1.size))
  sh2.orders <- matrix(NA,iters,nrow(DATA)-ncol(sh1.orders))
  sh.dets <- matrix(NA,iters,2)
  score.cors <- loadings.cors <- array(NA,dim=c(max.rank,max.rank,iters)) ## this is the maximum size it could be...


  pred.fi.array <- pred.u.array <- array(NA,dim=c(nrow(DATA),min(dim(DATA)),iters))
  pred.sds <- pred.mds <- matrix(NA,nrow(DATA),iters)
  rownames(pred.sds) <- rownames(pred.mds) <- rownames(DATA)


  for(i in 1:iters){

    ## the sort is simply for aesthetics...
    sh1 <- sort(sample(nrow(DATA),ceiling(nrow(DATA)*sh1.size)))
    sh2 <- sort(setdiff(1:nrow(DATA),sh1))

    sh1.data <- expo.scale(DATA[sh1,],center=center,scale=scale)
    sh1.res <- gsvd(sh1.data,k=min(k,max.rank))
    sh1.center <- attributes(sh1.data)$`scaled:center`
    sh1.scale <- attributes(sh1.data)$`scaled:scale`
    rm(sh1.data) # help the memory footprint

    sh2.data <- expo.scale(DATA[sh2,],center=center,scale=scale)
    sh2.res <- gsvd(sh2.data,k=min(k,max.rank))
    sh2.center <- attributes(sh2.data)$`scaled:center`
    sh2.scale <- attributes(sh2.data)$`scaled:scale`
    rm(sh2.data) # help the memory footprint


    sh.dets[i,] <- c(geometric.mean(sh1.res$d^2),geometric.mean(sh2.res$d^2))
    sh1.orders[i,] <- sh1
    sh2.orders[i,] <- sh2

    loadings.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <-
      cor(
        sh1.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))],
        sh2.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))]
      )


    ##### THIS BLOCK NEEDS FIXIN'

    ## predict based on the center/scale of the OTHER half.
      ### this should be fixed here, too...

        ## this block should only produce the eventual distances, we don't need the actual scores; just correlations and MD & SD.

    sh1.pred.fi <- expo.scale(DATA[sh1,],center=sh2.center,scale=sh2.scale) %*% sh2.res$v
    sh1.pred.u <- sweep(sh1.pred.fi,2,sh2.res$d,"/")

    sh2.pred.fi <- expo.scale(DATA[sh2,],center=sh1.center,scale=sh1.scale) %*% sh1.res$v
    sh2.pred.u <- sweep(sh2.pred.fi,2,sh1.res$d,"/")

    score.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <-
      (
        (cor(
          sh1.res$u[,1:min(ncol(sh1.res$u), ncol(sh1.pred.u))],
          sh1.pred.u[,1:min(ncol(sh1.res$u), ncol(sh1.pred.u))]
      )^2)[,1: min(ncol(sh1.res$u), ncol(sh1.pred.u),ncol(sh2.res$u), ncol(sh2.pred.u)) ] +
        (cor(
          sh2.res$u[,1:min(ncol(sh2.res$u), ncol(sh2.pred.u))],
          sh2.pred.u[,1:min(ncol(sh2.res$u), ncol(sh2.pred.u))]
      )^2)[,1: min(ncol(sh1.res$u), ncol(sh1.pred.u),ncol(sh2.res$u), ncol(sh2.pred.u)) ]
      ) /2



    pred.sds[sh1, i] <- rowSums(sh1.pred.fi^2)
    pred.sds[sh2, i] <- rowSums(sh2.pred.fi^2)

    pred.mds[sh1, i] <- rowSums(sh1.pred.u^2)
    pred.mds[sh2, i] <- rowSums(sh2.pred.u^2)

      ## these could be used for distances derived from subspaces...
      pred.fi.array[sh1, 1:min(c(ncol(pred.fi.array), ncol(sh1.pred.fi))) ,i] <- sh1.pred.fi
      pred.fi.array[sh2, 1:min(c(ncol(pred.fi.array), ncol(sh2.pred.fi))) ,i] <- sh2.pred.fi

      pred.u.array[sh1, 1:min(c(ncol(pred.u.array), ncol(sh1.pred.u))) ,i] <- sh1.pred.u
      pred.u.array[sh2, 1:min(c(ncol(pred.u.array), ncol(sh2.pred.u))) ,i] <- sh2.pred.u

    ##### THIS BLOCK NEEDS FIXIN'

    #print(i)
  }

  ## all those distances can be computed here.

  return( list(pred.sds=pred.sds,pred.mds=pred.mds, pred.fi.array=pred.fi.array,pred.u.array=pred.u.array,sh1.orders=sh1.orders,sh2.orders=sh2.orders,sh.dets=sh.dets,loadings.cors=loadings.cors,score.cors=score.cors) )
}

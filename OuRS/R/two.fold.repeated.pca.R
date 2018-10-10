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
  rownames(pred.fi.array) <- rownames(pred.u.array) <- rownames(DATA)


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

    ## predict based on the center/scale of the OTHER half.
    pred.fi.array[sh1,1:length(sh2.res$d),i] <- expo.scale(DATA[sh1,],center=sh2.center,scale=sh2.scale) %*% sh2.res$v
    pred.u.array[sh1,1:length(sh2.res$d),i] <- sweep(pred.fi.array[sh1,1:length(sh2.res$d),i],2,sh2.res$d,"/")

    pred.fi.array[sh2,1:length(sh1.res$d),i] <- expo.scale(DATA[sh2,],center=sh1.center,scale=sh1.scale) %*% sh1.res$v
    pred.u.array[sh2,1:length(sh1.res$d),i] <- sweep(pred.fi.array[sh2,1:length(sh1.res$d),i],2,sh1.res$d,"/")

    # the ODs can be brought back here by projecting onto each others subspaces.


    ## we can have so many bells and whistles...
    loadings.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <- cor(sh1.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))],sh2.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))])


    score.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <- cor(rbind(sh1.res$u[,1:min(c(length(sh1.res$d),length(sh2.res$d)))],sh2.res$u[,1:min(c(length(sh1.res$d),length(sh2.res$d)))]),pred.u.array[c(sh1,sh2),1:min(c(length(sh1.res$d),length(sh2.res$d))),i])

    #print(i)
  }

  ## all those distances can be computed here.

  return( list(pred.fi.array=pred.fi.array,pred.u.array=pred.u.array,sh1.orders=sh1.orders,sh2.orders=sh2.orders,sh.dets=sh.dets,loadings.cors=loadings.cors,score.cors=score.cors) )
}

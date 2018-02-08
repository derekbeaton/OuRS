split.half.pca <- function(DATA,center=T,scale=F,iters=500,sh1.k=0,sh2.k=0){

  sh1.orders <- matrix(NA,iters,ceiling(nrow(DATA)/2))
  sh2.orders <- matrix(NA,iters,nrow(DATA)-ncol(sh1.orders))
  sh.dets <- matrix(NA,iters,2)
  loadings.cors <- array(NA,dim=c(min(dim(DATA)),min(dim(DATA)),iters)) ## this is the maximum size it could be...


  pred.fi.array <- pred.u.array <- array(NA,dim=c(nrow(DATA),min(dim(DATA)),iters))
  rownames(pred.fi.array) <- rownames(pred.u.array) <- rownames(DATA)


  for(i in 1:iters){

    ## the sort is simply for aesthetics...
    sh1 <- sort(sample(nrow(DATA),ceiling(nrow(DATA)/2)))
    sh2 <- sort(setdiff(1:nrow(DATA),sh1))

    sh1.data <- expo.scale(DATA[sh1,],center=center,scale=scale)
    sh1.res <- gsvd(sh1.data,k=sh1.k)
    sh1.center <- attributes(sh1.data)$`scaled:center`
    sh1.scale <- attributes(sh1.data)$`scaled:scale`
    rm(sh1.data) # help the memory footprint

    sh2.data <- expo.scale(DATA[sh2,],center=center,scale=scale)
    sh2.res <- gsvd(sh2.data,k=sh2.k)
    sh2.center <- attributes(sh2.data)$`scaled:center`
    sh2.scale <- attributes(sh2.data)$`scaled:scale`
    rm(sh2.data) # help the memory footprint

    ## we can have so many bells and whistles...
    loadings.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <- cor(sh1.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))],sh2.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))])

    sh.dets[i,] <- c(geometric.mean(sh1.res$d^2),geometric.mean(sh2.res$d^2))
    sh1.orders[i,] <- sh1
    sh2.orders[i,] <- sh2

    ## predict based on the center/scale of the OTHER half.
    pred.fi.array[sh1,1:length(sh2.res$d),i] <- expo.scale(DATA[sh1,],center=sh2.center,scale=sh2.scale) %*% sh2.res$v
    pred.u.array[sh1,1:length(sh2.res$d),i] <- sweep(pred.fi.array[sh1,1:length(sh2.res$d),i],2,sh2.res$d,"/")

    pred.fi.array[sh2,1:length(sh1.res$d),i] <- expo.scale(DATA[sh2,],center=sh1.center,scale=sh1.scale) %*% sh1.res$v
    pred.u.array[sh2,1:length(sh1.res$d),i] <- sweep(pred.fi.array[sh2,1:length(sh1.res$d),i],2,sh1.res$d,"/")

    print(i)
  }


  return( list(pred.fi.array=pred.fi.array,pred.u.array=pred.u.array,sh1.orders=sh1.orders,sh2.orders=sh2.orders,sh.dets=sh.dets,loadings.cors=loadings.cors) )
}

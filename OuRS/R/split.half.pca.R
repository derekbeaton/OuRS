split.half.pca <- function(DATA,center=T,scale=F,iters=500,k=0){

  return(two.fold.repeated.pca(DATA=DATA,center=center,scale=scale,iters=iters,sh1.size=.5,k=k))

}

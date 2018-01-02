  ## needs a better name.
quant.mcd.sim <- function(data,best.sample,iters=100){

  for(i in 1:iters){
    mvrnorm(nrow(data),colMeans(data[best.sample,]),cov(data[best.sample,]))
  }

}

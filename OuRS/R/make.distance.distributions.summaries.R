make.distance.distributions.summaries <- function(distance.array,components=1:dim(distance.array)[2],types="all",lower.percentile=.025,upper.percentile=.975){

  distance.array[is.na(distance.array)] <- 0

      ## allows us to compute partial if we'd like.
  these.dists <- sqrt(apply(distance.array[,components,]^2,c(1,3),sum))


    ## I ain't got time for you to not read the documentation!
  if(length(types)!=1){
    types <- "all"
  }
  if(!(types %in% c("all","median","iqr","percentile"))){
    types <- "all"
  }

  compute.percentiles <- compute.iqr <- compute.median <- F
  if(types==tolower("all")){
    compute.percentiles <- compute.iqr <- compute.median <- T
  }
  if(types==tolower("median")){
    compute.median <- T
  }
  if(types==tolower("iqr")){
    compute.iqr <- T
  }
  if(types==tolower("percentile")){
    compute.median <- T
  }



  if(compute.median){
    median.dist <- apply(these.dists,1,median)
  }else{
    median.dist <- NULL
  }
  if(compute.iqr){
    iqr.dist <- apply(these.dists,1,IQR)
  }else{
    iqr.dist <- NULL
  }
  if(compute.percentiles){
  percentile.dist <- apply(these.dists,1,
        function(x){
          sort(x)[ceiling(length(x)*upper.percentile)] - sort(x)[floor(length(x)*lower.percentile)]
        })
  }else{
    percentile.dist <- NULL
  }

  return(list(dists=these.dists,median.dist=median.dist,iqr.dist=iqr.dist,percentile.dist=percentile.dist))

}

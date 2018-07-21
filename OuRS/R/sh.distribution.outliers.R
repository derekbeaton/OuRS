sh.distribution.outliers <- function(distances, alpha=.75){

  outlier.threshold <- sort(c(dist.distrs))[round( length(c(dist.distrs)) * alpha )]

  outlier.scores <- apply(distances,1,function(x){sum(x >= upper.bound) / length(x)})

  return(list(outlier.threshold=outlier.threshold, outlier.scores=outlier.scores))

}

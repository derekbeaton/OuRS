sh.distribution.outliers <- function(distances, alpha=.75){

  outlier.threshold <- sort(c(distances))[round( length(c(distances)) * alpha )]

  outlier.scores <- apply(distances,1,function(x){sum(x >= outlier.threshold) / length(x)})

  return(list(outlier.threshold=outlier.threshold, outlier.scores=outlier.scores))

}

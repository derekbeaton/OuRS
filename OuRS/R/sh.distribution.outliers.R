sh.distribution.outliers <- function(distances, distribution.alpha=.75, outlier.alpha = distribution.alpha){

  outlier.threshold <- sort(c(distances))[round( length(c(distances)) * distribution.alpha )]
  outlier.scores <- apply(distances,1,function(x){sum(x >= outlier.threshold) / length(x)})

  outliers <- outlier.scores >= outlier.alpha

  return(list(outlier.threshold=outlier.threshold, outlier.scores=outlier.scores, outliers = outliers))

}

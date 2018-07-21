dist.array.outliers <- function(dist.array, total.dist.cutoff = .95, outlier.cutoff = .95){

  dist.array[is.na(dist.array)] <- 0

  dist.distrs <- sqrt(apply(dist.array^2,c(1,3),sum))

  upper.bound <- sort(c(dist.distrs))[round( length(c(dist.distrs)) * total.dist.cutoff )]
  outlier.scores <- apply(dist.distrs,1,function(x){sum(x >= upper.bound) / length(x)})
  outlier.threshold <- outlier.scores > outlier.cutoff

  return(list(dists=dist.distrs, cut.point = upper.bound, outlier.probabilities = outlier.scores, outliers = outlier.threshold))

}

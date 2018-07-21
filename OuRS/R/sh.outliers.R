sh.outliers <- function(sh.out, total.dist.cutoff = .95, outlier.cutoff = .95){

  score.outliers <- dist.array.outliers(sh.out$pred.fi.array, total.dist.cutoff = total.dist.cutoff, outlier.cutoff = outlier.cutoff)
  mahal.outliers <- dist.array.outliers(sh.out$pred.u.array, total.dist.cutoff = total.dist.cutoff, outlier.cutoff = outlier.cutoff)

  return(list(score.outliers=score.outliers,mahal.outliers=mahal.outliers))

}

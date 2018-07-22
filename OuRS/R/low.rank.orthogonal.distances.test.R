low.rank.orthogonal.distances.test <- function(DATA, center=T, scale=F, components=1:2, bootstrap.iters = 100, alpha = .75, bootstrap.shortcut = T){

  DATA <- expo.scale(DATA,center=center,scale=scale)
  full.svd.res <- gsvd(DATA)

  if(length(components>1)){
    low.rank.rebuild <- full.svd.res$u[,components] %*% diag(full.svd.res$d[components]) %*% t(full.svd.res$v[,components])
  }else if(length(components)==1){
    low.rank.rebuild <- as.matrix(full.svd.res$u[,components]) %*% as.matrix(full.svd.res$d[components]) %*% t(as.matrix(full.svd.res$v[,components]))
  }else{
    warning("Invalid number of components selected. First two will be used.")
    components <- 1:2
    low.rank.rebuild <- full.svd.res$u[,components] %*% diag(full.svd.res$d[components]) %*% t(full.svd.res$v[,components])
  }

  od <- sqrt(rowSums( (DATA - low.rank.rebuild)^2 ))

  if(bootstrap.iters > 0){

    if(bootstrap.shortcut){

      all.ods <- sample(od,size=length(od)*bootstrap.iters,replace=T)

    }else{

      all.ods.mat <- matrix(NA,length(od),bootstrap.iters)
      for(i in 1:bootstrap.iters){
        boot.samp <- sample(length(od), length(od), replace=T)
        boot.DATA <- expo.scale(DATA[boot.samp,],center=attributes(DATA)$`scaled:center`,scale=attributes(DATA)$`scaled:scale`)


          if(length(components>1)){
            boot.low.rank.rebuild <- (boot.DATA %*% full.svd.res$v[,components]) %*% t(full.svd.res$v[,components])
          }else {
            boot.low.rank.rebuild <- (boot.DATA %*% as.matrix(full.svd.res$v[,components])) %*% t(as.matrix(full.svd.res$v[,components]))
          }

        all.ods.mat[,i] <- sqrt(rowSums( (boot.DATA - boot.low.rank.rebuild)^2 ))

      }
      all.ods <- c(all.ods.mat)
      rm(all.ods.mat)


    }
    outlier.threshold <- sort(all.ods)[ceiling(length(all.ods)*alpha)]
    outliers <- od >= outlier.threshold
    return(list(od=od, outliers=outliers, outlier.threshold=outlier.threshold))
  }

  return(list(od=od, outliers=rep(FALSE,length(od)), outlier.threshold=0))


}

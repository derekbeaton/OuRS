low.rank.orthogonal.distances <- function(DATA, center=T, scale=F, components=1:2){

  DATA <- expo.scale(DATA,center=center,scale=scale)
  full.svd.res <- gsvd(DATA)

  if(length(components>1)){
    low.rank.rebuild <- full.svd.res$u[,components] %*% diag(full.svd.res$d[components]) %*% t(full.svd.res$v[,components])
  }else if(length(components)==1){
    low.rank.rebuild <- as.matrix(full.svd.res$u[,components]) %*% as.matrix(full.svd.res$d[components]) %*% t(as.matrix(full.svd.res$v[,components]))
  }else{
    warning("Invalid number of components selected. First two will be used.")
    components <- 1:2
  }
  S <- DATA - low.rank.rebuild
  od <- sqrt(rowSums(s.mat^2))

  return(list(S=S,od=od))

}

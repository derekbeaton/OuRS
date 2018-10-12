
## update to CA. stolen primarily from SlimPosition (https://github.com/derekbeaton/ExPosition-Family/blob/master/ExPosition2/SlimPosition/Package/R/sp.ca.R)

ca <- function(DATA, k = 0, compact = F){

  sum.data <- sum(DATA)
  rowSums.data <- rowSums(DATA)
  wi <- rowSums.data/sum.data
  wj <- colSums(DATA)/sum.data

  res <- gsvd( sweep(sweep(DATA,1,rowSums.data,"/"),2,wj), wi, 1/wj, k = k )
  res$fi <- sweep(res$fi,1,wi,"/")
  if(compact){
    res <- list(fi=res$fi, fj=res$fj, d.orig=res$d.orig, u=res$u, v=res$v)
  }
  res$type <- "ca"
  return(res)
}

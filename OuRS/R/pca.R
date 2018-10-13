## Stolen primarily from SlimPosition (https://github.com/derekbeaton/ExPosition-Family/blob/master/ExPosition2/SlimPosition/Package/R/sp.pca.R)
## Compact will be only: fi, fj, u, v, and d.orig
pca <- function(DATA, center = T, scale = T, k = 0, compact = F){

  res <- gsvd(expo.scale(DATA, center = center, scale = scale), k = k)
  if(compact){
    res <- list(fi=res$fi, fj=res$fj, d.orig=res$d.orig, u=res$u, v=res$v)
  }
  res$type <- "pca"
  return(res)
}
